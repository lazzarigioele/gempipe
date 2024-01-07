import os
import pickle
import multiprocessing
import itertools
import shutil
import subprocess


import pandas as pnd
import cobra
from cobra.util.solver import linear_reaction_coefficients
from Bio import SeqIO, SeqRecord, Seq


from ..commons import get_blast_header
from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_retained_accessions



def task_brh(proteome, args):
    
    
    # retrive the arguments:
    ref_proteome = args['ref_proteome']
    
    
    # get the basename without extension:
    basename = os.path.basename(proteome)
    accession, _ = os.path.splitext(basename)
    
    
    # create subdir without overwriting: 
    os.makedirs(f'working/brh/{accession}/', exist_ok=True)
    os.makedirs(f'working/brh/{accession}/dbs/reference/', exist_ok=True)
    os.makedirs(f'working/brh/{accession}/dbs/{accession}/', exist_ok=True)
            
    
    # create blast database for reference: 
    shutil.copyfile(ref_proteome, f'working/brh/{accession}/dbs/reference/ref_proteome.faa')  # just the content, not the permissions.  
    command = f"""makeblastdb -in working/brh/{accession}/dbs/reference/ref_proteome.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()

    
    # create blast database for current strain:
    shutil.copyfile(proteome, f'working/brh/{accession}/dbs/{accession}/{accession}.faa')  # just the content, not the permissions.    
    command = f"""makeblastdb -in working/brh/{accession}/dbs/{accession}/{accession}.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    

    # perform blastp for reference-on-accession: 
    command = f'''blastp \
        -query working/brh/{accession}/dbs/reference/ref_proteome.faa \
        -db working/brh/{accession}/dbs/{accession}/{accession}.faa \
        -out working/brh/{accession}/align_ref_vs_acc.tsv \
        -outfmt "6 {get_blast_header()}"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # perform blastp for accession-on-reference: 
    command = f'''blastp \
        -query working/brh/{accession}/dbs/{accession}/{accession}.faa \
        -db working/brh/{accession}/dbs/reference/ref_proteome.faa \
        -out working/brh/{accession}/align_acc_vs_ref.tsv \
        -outfmt "6 {get_blast_header()}"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # read both the alignments:
    header = f"{get_blast_header()}".split(' ')
    align_ref_vs_acc = pnd.read_csv(f'working/brh/{accession}/align_ref_vs_acc.tsv', names=header, sep='\t')
    align_ref_vs_acc['qcov'] = round((align_ref_vs_acc['qend'] -  align_ref_vs_acc['qstart'] +1)/ align_ref_vs_acc['qlen'] * 100, 1)
    align_acc_vs_ref = pnd.read_csv(f'working/brh/{accession}/align_acc_vs_ref.tsv', names=header, sep='\t')
    align_acc_vs_ref['qcov'] = round((align_acc_vs_ref['qend'] -  align_acc_vs_ref['qstart'] +1)/ align_acc_vs_ref['qlen'] * 100, 1)
    

    # parse the alignments
    results_df = [] 
    for cds in align_acc_vs_ref['qseqid'].unique():

        
        # acc_vs_ref
        curr_hsps = align_acc_vs_ref[align_acc_vs_ref['qseqid'] == cds]
        curr_hsps_filt = curr_hsps[(curr_hsps['qcov'] >= 70) & (curr_hsps['evalue'] <= 1e-5)]
        curr_hsps_filt_sort = curr_hsps_filt.sort_values(by='evalue', ascending=True)
        curr_hsps_filt_sort.reset_index(inplace=True, drop=True)
        try: best_hit = curr_hsps_filt_sort.loc[0, 'sseqid']
        except: 
            results_df.append({'cds': cds, 'ref': None, 'reciprocal': 'NA'})
            continue


        # ref_vs_acc
        curr_hsps2 = align_ref_vs_acc[align_ref_vs_acc['qseqid'] == best_hit]
        curr_hsps_filt2 = curr_hsps2[(curr_hsps2['qcov'] >= 70) & (curr_hsps2['evalue'] <= 1e-5)]
        curr_hsps_filt_sort2 = curr_hsps_filt2.sort_values(by='evalue', ascending=True)
        curr_hsps_filt_sort2.reset_index(inplace=True, drop=True)
        try: best_hit2 = curr_hsps_filt_sort2.loc[0, 'sseqid']
        except: 
            results_df.append({'cds': cds, 'ref': best_hit, 'reciprocal': 'NA'})
            continue

        
        # annotate if bi-directional or mono-directional
        if cds == best_hit2: results_df.append({'cds': cds, 'ref': best_hit, 'reciprocal': '<=>'})
        else: results_df.append({'cds': cds, 'ref': best_hit, 'reciprocal': '=>'})


    # save results to disk
    ref_proteome_basename = os.path.basename(ref_proteome)
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv(f'working/brh/{accession}_brh_{ref_proteome_basename}.csv')
    
    
    # save disk space removeing databases: 
    shutil.rmtree(f'working/brh/{accession}/') 
    
    
    # return a row for the dataframe
    return [{'accession': accession, 'completed': True}]
    


def create_refgid_to_clusters(refmodel_basename, ref_proteome_basename): 
    
    
    # load the previously created doctionaries: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    with open('working/clustering/seq_to_cluster.pickle', 'rb') as handler:
        seq_to_cluster = pickle.load(handler)
    
    
    # parse each filtered accession: 
    refgid_to_clusters = {}
    for species in species_to_proteome.keys():
        for proteome in species_to_proteome[species]: 
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
    
            
            # read the brh results for this accession: 
            df_result = pnd.read_csv(f'working/brh/{accession}_brh_{ref_proteome_basename}.csv', index_col=0)
            df_result = df_result.set_index('cds', drop=True, verify_integrity=True)
            
            
            # get only the reciprocal hits: 
            df_brh = df_result[df_result['reciprocal'] == '<=>']
            
            
            # populate the dictionary: 
            for cds, row in df_brh.iterrows(): 
                if row['ref'] not in refgid_to_clusters.keys(): 
                    refgid_to_clusters[row['ref']] = set()
                refgid_to_clusters[row['ref']].add(seq_to_cluster[cds])


    # save the dictionary: 
    with open(f'working/brh/{refmodel_basename}.refgid_to_clusters.pickle', 'wb') as handler:
        pickle.dump(refgid_to_clusters, handler)
        
        
        
def translate_refmodel(logger, refmodel, ref_proteome): 
    
    
    # set up the cobra solver
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"

    
    # load the model according to the file type
    logger.info("Loading the provided reference model...")
    refmodel_basename = os.path.basename(refmodel)
    if refmodel.endswith('.json'): 
        refmodel = cobra.io.load_json_model(refmodel)
    elif refmodel.endswith('.sbml'): 
        refmodel = cobra.io.read_sbml_model(refmodel)
    elif refmodel.endswith('.xml'): 
        refmodel = cobra.io.read_sbml_model(refmodel)
    else:
        logger.error(refmodel + ": extension not recognized.")
        return 1 
    logger.info(f"Done, {' '.join(['G:', str(len(refmodel.genes)), '|', 'R:', str(len(refmodel.reactions)), '|', 'M:', str(len(refmodel.metabolites))])}.")
    
    
    # print the preloaded objective reaction
    objs = list(linear_reaction_coefficients(refmodel).keys())
    if len(objs) > 1: logger.warning("More than 1 objective reactions were set up. Showing the first.")
    logger.info(f"The following objective was set up: {objs[0].id}.")
    
    
    # save the reference model in a standard format
    logger.debug("Saving a copy of the reference model in JSON format...")
    cobra.io.save_json_model(refmodel, f'working/brh/{refmodel_basename}.refmodel_original.json') # ext can be repeated.
    
    
    # create a copy, later translated
    logger.info("Converting reference model's gene notation to clusters...")
    refmodel_t = refmodel.copy()

    
    # get the modeled genes: 
    modeled_gids = set([g.id for g in refmodel.genes])


    # sometimes the reference proteome contains less genes respect to those modeled.
    # we remove the reference genes missing from the proteome: 
    gids_in_proteome = set()
    with open(ref_proteome, 'r') as r_handler:                  
        for seqrecord in SeqIO.parse(r_handler, "fasta"):
            gid = seqrecord.id
            gids_in_proteome.add(gid)
    # remove non-modeled genes: 
    if len(modeled_gids - gids_in_proteome) > 0:
        to_remove = list(modeled_gids - gids_in_proteome)
        logger.info(f"The following genes will be removed from the reference model, as they do not appear in the reference proteome: {list(modeled_gids - gids_in_proteome)}.") 
        cobra.manipulation.remove_genes(refmodel_t, to_remove, remove_reactions=True)


    # load the refgid_to_clusters dictionary (1-to-many)
    with open(f'working/brh/{refmodel_basename}.refgid_to_clusters.pickle', 'rb') as handler:
        refgid_to_clusters = pickle.load(handler)

                
    # finally rename the genes: 
    for r in refmodel_t.reactions:
        if any([gid in r.gene_reaction_rule for gid in gids_in_proteome]): 
            gpr = r.gene_reaction_rule
            # force each gid to be surrounded by spaces: 
            gpr = ' ' + gpr.replace('(', ' ( ').replace(')', ' ) ') + ' '

            
            # translate this GPR
            for gid in refgid_to_clusters.keys():
                if f' {gid} ' in gpr:  
                    gpr = gpr.replace(f' {gid} ', f' ({" or ".join(refgid_to_clusters[gid])}) ')

                    
            # remove spaces between parenthesis
            gpr = gpr.replace(' ( ', '(').replace(' ) ', ')')
            # remove spaces at the extremes: 
            gpr = gpr[1: -1]
            
            
            # New genes are introduced. Parethesis at the extremes are removed if not necessary. 
            r.gene_reaction_rule = gpr
            r.update_genes_from_gpr()


    # now remove the reference genes:
    to_remove = [g for g in refmodel_t.genes if g.id in gids_in_proteome]
    cobra.manipulation.delete.remove_genes(refmodel_t, to_remove, remove_reactions=True)
    logger.info(f"Done, {' '.join(['G:', str(len(refmodel_t.genes)), '|', 'R:', str(len(refmodel_t.reactions)), '|', 'M:', str(len(refmodel_t.metabolites))])}.")
    
    
    # save the reference model in a standard format
    logger.debug("Saving a copy of the converted reference model in JSON format...")
    cobra.io.save_json_model(refmodel_t, f'working/brh/{refmodel_basename}.refmodel_translated.json')  # ext can be repeated.
    
    
    return 0
    
    

def perform_brh(logger, cores, ref_proteome): 
    
    
    # some log messages:
    logger.info("Performing the best reciprocal hits (BRH) alignment against the reference proteome...")
    if not os.path.exists(ref_proteome): # check the input:
        logger.error(f"Provided path to the reference proteome (-rp/--ref_proteome) does not exist: {ref_proteome}.")
        return 1
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/brh/', exist_ok=True)

    
    # load the previously created species_to_genome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    

    # check if it's everything pre-computed
    ref_proteome_basename = os.path.basename(ref_proteome)
    results_presence = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            results_presence.append(os.path.exists(f'working/brh/{accession}_brh_{ref_proteome_basename}.csv'))
    if all(results_presence): 
        # log some message: 
        logger.info('Found all the needed files already computed. Skipping this step.')
        # signal to skip this module:
        return 0
        
    

    # create items for parallelization: 
    items = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]: 
            items.append(proteome)
    

    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession', 'completed']), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_brh), 
            itertools.repeat({'ref_proteome': ref_proteome}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    return 0



def convert_reference(logger, refmodel, ref_proteome):
    
    
    # some log messages:
    logger.info("Translating the reference model's genes to clusters...")
    if not os.path.exists(refmodel): # check the input:
        logger.error(f"Provided path to the reference model (-rm/--ref_model) does not exist: {refmodel}.")
        return 1
    
    
    # get basename for reference model:
    refmodel_basename = os.path.basename(refmodel)
    
    
    # check if it's everything pre-computed
    if os.path.exists('working/brh/proc_acc.pickle'):
        with open('working/brh/proc_acc.pickle', 'rb') as handler:
            proc_acc = pickle.load(handler) 
        if get_retained_accessions() == proc_acc:
            if os.path.exists(f'working/brh/{refmodel_basename}.refmodel_original.json'):
                if os.path.exists(f'working/brh/{refmodel_basename}.refmodel_translated.json'):
                    if os.path.exists(f'working/brh/{refmodel_basename}.refgid_to_clusters.pickle'):
                        # log some message: 
                        logger.info('Found all the needed files already computed. Skipping this step.')
                        # signal to skip this module:
                        return 0

    
    # create a dictionary ref_seq-to-clusters, parsing the BRHs. 
    ref_proteome_basename = os.path.basename(ref_proteome)
    create_refgid_to_clusters(refmodel_basename, ref_proteome_basename)
    
    
    # get a opy of the refmodel, and translate its genes to clusters notation. 
    response = translate_refmodel(logger, refmodel, ref_proteome)
    if response == 1: return 1


    # make traces to keep track of the accessions processed:
    shutil.copyfile('working/annotation/proc_acc.pickle', 'working/brh/proc_acc.pickle')
    
    
    return 0
    