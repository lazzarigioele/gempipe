import multiprocessing
import itertools
import os
import shutil


import cobra


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import fba_no_warnings



def subtract_clusters(ss_model, panmodel, pam, accession):
    
    
    modeled_gids = [g.id for g in panmodel.genes]  # get medoled genes ID
    to_remove = []  # genes to delete
    
    
    # iterate the PAM :
    for cluster in pam.index: 
        # consider only if it is modeled:
        if cluster in modeled_gids: 
            cell = pam.loc[cluster, accession]
            if type(cell) == float:  # empty cell
                to_remove.append(ss_model.genes.get_by_id(cluster))
                continue
            # get all the sequences not containing a premature stop:
            seqs = [i for i in cell.split(';') if i.endswith('_stop')==False]
            if len(seqs) == 0:  # they were all '_stop' sequences.
                to_remove.append(ss_model.genes.get_by_id(cluster))
                continue
    
    
    # delete marked genes
    cobra.manipulation.delete.remove_genes(ss_model, to_remove, remove_reactions=True)
    
    
    return 0



def translate_genes(logger, ss_model, panmodel, pam, accession):
    
    
    # create a dups dict: 
    dups_dict = {}  # {'Cluster_XXXX': ['GFDSA_XXXX', 'GFDSA_XXXX', 'GFDSA_XXXX']}
    modeled_gids = [g.id for g in panmodel.genes]  # get medoled genes ID
    
    
    # iterate the PAM :
    for cluster in pam.index: 
        # consider only if it is modeled:
        if cluster in modeled_gids: 
            cell = pam.loc[cluster, accession]
            if type(cell) != float:  
                # get all the sequences not containing a premature stop:
                seqs = [i for i in cell.split(';') if i.endswith('_stop')==False]
                if len(seqs) >= 1:  
                    # populate the dict:
                    dups_dict[cluster] = seqs
    
    
    # iterate reactions:
    for r in ss_model.reactions:
        if 'Cluster_' in r.gene_reaction_rule: 
            gpr = r.gene_reaction_rule
            
            
            # force each gid to be surrounded by spaces: 
            gpr = ' ' + gpr.replace('(', ' ( ').replace(')', ' ) ') + ' '
            for cluster in dups_dict.keys():
                
                
                # search this gid surrounded by spaces:
                if f' {cluster} ' in gpr:
                    gpr = gpr.replace(f' {cluster} ', f' ({" or ".join(dups_dict[cluster])}) ')

            
            # remove spaces between parenthesis
            gpr = gpr.replace(' ( ', '(').replace(' ) ', ')')
            # remove spaces at the extremes: 
            gpr = gpr[1: -1]
            
            
            # New genes are introduced. Parethesis at the extremes are removed if not necessary. 
            r.gene_reaction_rule = gpr
            r.update_genes_from_gpr()
            
            
    # remaining old 'Cluster_'s need to removed.
    # remove if (1) hte ID starts with clusters AND (2) they are no more associated with any reaction
    to_remove = [g for g in ss_model.genes if (g.id.startswith("Cluster_") and len(g.reactions)==0)]
    cobra.manipulation.delete.remove_genes(ss_model, to_remove, remove_reactions=True)
    
    
    # check if the cleaning was successful
    for g in ss_model.genes:
        if g.id.startswith("Cluster_"):
            logger.error(f"Gene ID translation was not successful ({g.id} remained, accession {accession}). Please contact the developer.")
            return 1
        
        
    return 0



def annotate_genes(ss_model, accession, gannots):
    
    
    
    if gannots is not None: 
        # select the rows belonging to this accession: 
        gannots_accession = gannots[gannots['accession']==accession]
        # get the modeled gids: 
        gids_modeled = set([g.id for g in ss_model.genes])


        # iterate the 'gannots' dataframe:
        for index, row in gannots_accession.iterrows(): 
            if row['gid'] in gids_modeled:
                g = ss_model.genes.get_by_id(row['gid'])


                # '== str' to exclude empty cells: 
                if type(row['refseq']) == str:
                    g.annotation['refseq'] = [row['refseq']]
                if type(row['ncbiprotein']) == str:
                    g.annotation['ncbiprotein'] = [row['ncbiprotein']]
                if type(row['ncbigene']) == str:
                    g.annotation['ncbigene'] = [row['ncbigene']]
                if type(row['kegg']) == str:
                    g.annotation['kegg.genes'] = [row['kegg']]
                if type(row['uniprot']) == str:
                    g.annotation['uniprot'] = [row['uniprot']]
                
                
    # SBO annotation has been lost during the gene ID translation:
    for g in ss_model.genes:
        g.annotation['sbo'] = ['SBO:0000243']  # generic gene
                

def task_derivestrain(accession, args):
    
    
    # get the arguments
    panmodel = args['panmodel']
    pam = args['pam']
    report = args['report']
    outdir = args['outdir']
    gannots = args['gannots']
    logger = args['logger']
    sbml = args['sbml']
    
    
    # define key objects: 
    ss_model = panmodel.copy()  # create strain specific model
    
    
    # remove Clusters
    response = subtract_clusters(ss_model, panmodel, pam, accession)
    if response == 1: return 1   # will raise an error since a dict is expected
    
    # translate gene IDs:
    response = translate_genes(logger, ss_model, panmodel, pam, accession)
    if response == 1: return 1   # will raise an error since a dict is expected
    
    # annotate genes:
    response = annotate_genes(ss_model, accession, gannots)
    if response == 1: return 1   # will raise an error since a dict is expected
    
    
    # get the associated species:
    report = report[report['accession'] == accession]
    if len(report) > 1: 
        logger.error("Duplicated accessions in the provided report. Please report this error to the developer.")
    report = report.reset_index(drop=True)
    species = report.loc[0, 'species']
    strain = report.loc[0, 'strain']
    niche = report.loc[0, 'niche']
    
    
    # get some metrics: 
    n_G = len(ss_model.genes)
    n_R = len(ss_model.reactions)
    n_M = len(ss_model.metabolites)
    
    
    # try the FBA: 
    res, obj_value, status = fba_no_warnings(ss_model)

    
    
    # save strain specific model to disk
    cobra.io.save_json_model(ss_model, f'{outdir}/{accession}.json')
    if sbml: cobra.io.write_sbml_model(ss_model, f'{outdir}/{accession}.xml')
    
    
    # compose the new row:
    return [{'accession': accession, 'species': species, 'strain': strain, 'niche': niche, 'G': n_G, 'R': n_R, 'M': n_M, 'obj_value': obj_value, 'status': status }]



def  derive_strain_specific(logger, outdir, cores, panmodel, pam, report, gannots, sbml):

    
    # log some messages
    logger.info("Deriving strain-specific models...")
    
   
    # create output dir
    if os.path.exists(outdir + 'strain_models/'):
        # always overwriting if already existing
        shutil.rmtree(outdir + 'strain_models/')  
    os.makedirs(outdir + 'strain_models/', exist_ok=True)
    

    # create items for parallelization: 
    items = []
    for accession in pam.columns:
        items.append(accession)
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession', 'species', 'strain', 'niche', 'G', 'R', 'M', 'obj_value', 'status']), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_derivestrain),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'panmodel': panmodel, 'pam': pam, 'report': report, 'outdir': outdir + 'strain_models', 'gannots': gannots, 'logger': logger, 'sbml': sbml}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save tabular output:
    df_derive_strains = all_df_combined
    df_derive_strains = df_derive_strains.sort_values(by=['species', 'strain', 'niche'], ascending=True)   # sort by species
    df_derive_strains.to_csv(outdir + 'derive_strains.csv')
    
    
    return 0