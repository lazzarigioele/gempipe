import os
import pickle
import multiprocessing
import itertools
import shutil
import subprocess
import warnings


import pandas as pnd
import cobra
from cobra.util.solver import linear_reaction_coefficients
from Bio import SeqIO, SeqRecord, Seq


from ..commons import get_blast_header
from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_retained_accessions
from ..commons import read_refmodel



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
    
    
    # save the alignments in a more readible formt:
    align_ref_vs_acc.to_csv(f'working/brh/{accession}_align_ref_vs_acc.csv')
    align_acc_vs_ref.to_csv(f'working/brh/{accession}_align_acc_vs_ref.csv')
    

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
    


def create_refgid_to_clusters(logger, refmodel_basename, ref_proteome_basename, edits_dict): 
    
    
    # load the previously created doctionaries: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    with open('working/clustering/seq_to_cluster.pickle', 'rb') as handler:
        seq_to_cluster = pickle.load(handler)
        
        
    # create a 'cds_to_newcluster' dict, parsing the edits_dict.
    # This is because one of the ss proteins in the brh results coul be 'broken'.
    # Being broken, its cluster would be changed.
    cds_to_newcluster = {}
    if edits_dict != None:
        # 'edits_dict' syntax is "x:cluster --> frag_x_y:newcluster"
        for key, value in edits_dict.items():
            cds = key.split(':', 1)[0]
            new_cluster = value.split(':', 1)[1]
            if cds in cds_to_newcluster.keys():
                logger.error("Each 'cds' should appear just one time in 'create_refgid_to_clusters()'.")
            cds_to_newcluster[cds] = new_cluster    
    # save the dictionary: 
    with open(f'working/brh/cds_to_newcluster.pickle', 'wb') as handler:
        pickle.dump(cds_to_newcluster, handler)
    
    
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
                    
                    
                # take into account eventual 'broken' proteins
                cluster = seq_to_cluster[cds]
                if edits_dict != None:
                    if cds in cds_to_newcluster.keys():
                        cluster = cds_to_newcluster[cds]
                    
                # finally add the cluster
                refgid_to_clusters[row['ref']].add(cluster)


    # save the dictionary: 
    with open(f'working/brh/{refmodel_basename}.refgid_to_clusters.pickle', 'wb') as handler:
        pickle.dump(refgid_to_clusters, handler)
        
        
        
def translate_refmodel(logger, refmodel, ref_proteome, refspont): 
    
    
    # load the model according to the file type
    logger.info("Loading the provided reference model...")
    refmodel_basename = os.path.basename(refmodel)
    refmodel = read_refmodel(refmodel)   # handle variuos formats
    if type(refmodel)==int: return 1  # an error was raised
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
    # we remove from the reference model the genes missing from the proteome: 
    gids_in_proteome = set()
    with open(ref_proteome, 'r') as r_handler:                  
        for seqrecord in SeqIO.parse(r_handler, "fasta"):
            gid = seqrecord.id
            gids_in_proteome.add(gid)
    # remove non-proteome genes (exclude the "spontaneous" gene): 
    if len(modeled_gids - gids_in_proteome) > 0:
        to_remove = []
        for i in list(modeled_gids - gids_in_proteome):
            if i != refspont:  # exclude the "spontaneous" gene
                to_remove.append(i)
        logger.info(f"The following genes will be removed from the reference model, as they do not appear in the reference proteome: {to_remove}.") 
        to_remove = [refmodel_t.genes.get_by_id(gid) for gid in to_remove]
        cobra.manipulation.remove_genes(refmodel_t, to_remove, remove_reactions=False)
        # update 'modeled_gids' set:
        modeled_gids = modeled_gids - set([g.id for g in to_remove])


    # load the refgid_to_clusters dictionary (1-to-many)
    with open(f'working/brh/{refmodel_basename}.refgid_to_clusters.pickle', 'rb') as handler:
        refgid_to_clusters = pickle.load(handler)
    
        
    # if gempipe is run in 'reference' mode, on a gempipe-generated reference pan-model,
    # the the old gene IDs and current clusteres will be conflicting. In these cases,
    # we add a prefix to the old gene ID, and correct 'modeled_gids' and 'refgid_to_clusters' accordingly.
    gempipe_on_gempipe = False
    if any(['Cluster_' in g.id for g in refmodel_t.genes]):
        gempipe_on_gempipe = True  # flag
        renaming_dict = {}
        for g in refmodel_t.genes:
            if g.id == refspont:
                continue
                
            # edit the 'modeled_gids' set:
            modeled_gids = modeled_gids - set([g.id])
            modeled_gids.add(f'old_{g.id}')
            
            # edit the 'gids_in_proteome' set:
            gids_in_proteome = gids_in_proteome - set([g.id])
            gids_in_proteome.add(f'old_{g.id}')
            
            # edit the 'refgid_to_clusters' dict:
            if g.id in refgid_to_clusters.keys():
                pointed = refgid_to_clusters[g.id]
                del refgid_to_clusters[g.id]
                refgid_to_clusters[f'old_{g.id}'] = pointed
                
            # populate the 'renaming_dict' (old: new)
            renaming_dict[g.id] = f'old_{g.id}'
            
        # finally rename genes with their prefixed version:
        cobra.manipulation.rename_genes(refmodel_t, renaming_dict)
            
        
    # at this point, the 'refgid_to_clusters' dict may not contain some reference genes.
    # This is because of the threshold applied earlier, on the 'qcov' and 'evalue', for the determination of the brh.
    # Here we try to recovere these genes: 
    translable_gids = set(list(refgid_to_clusters.keys()))
    no_brh_gids = (modeled_gids - translable_gids) - set([refspont])
    if len(no_brh_gids) > 0:
        logger.info(f"The genes {no_brh_gids} remained without a BRH. Now starting their recovery.")
        recovered_translations = {}
        
        # load the translation dict (also take into account the 'broken' proteins)
        with open('working/clustering/seq_to_cluster.pickle', 'rb') as handler:
            seq_to_cluster = pickle.load(handler)
        with open(f'working/brh/cds_to_newcluster.pickle', 'rb') as handler:
            cds_to_newcluster = pickle.load(handler)
        
        # iterate the 'good' accessions
        with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
            species_to_proteome = pickle.load(handler)
        for species in species_to_proteome.keys():
            for proteome in species_to_proteome[species]: 
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)

                # read the brh results for this accession: 
                ref_vs_acc = pnd.read_csv(f'working/brh/{accession}_align_ref_vs_acc.csv', index_col=0)
                
                # get the best match for each 'gid':
                for gid in no_brh_gids:
                    if gempipe_on_gempipe: gid = gid[len('old_'):]
                    matches = ref_vs_acc[ref_vs_acc['qseqid']==gid]
                    matches = matches.sort_values(by='evalue', ascending=True)
                    matches = matches.reset_index(drop=True)
                    
                    try: match = matches.iloc[0]
                    except: continue  # no match
                    
                    # populate the dictionary:
                    if gempipe_on_gempipe: gid = f'old_{gid}'
                    if gid not in recovered_translations.keys():
                        recovered_translations[gid] = set()
                    cds = match['sseqid']
                    cluster = seq_to_cluster[cds]
                    if cds in cds_to_newcluster.keys():  # manage also the case of 'broken' proteins
                        cluster = cds_to_newcluster[cds]
                    recovered_translations[gid].add(cluster) 
                    
        logger.debug(f'Recovered translations: {recovered_translations}')
        # supplement the 'refgid_to_clusters' dict:
        for key, value in recovered_translations.items():
            refgid_to_clusters[key] = value
                    
                
    # finally rename the genes: 
    for r in refmodel_t.reactions:
        gpr = r.gene_reaction_rule
        
        
        # force each gid to be surrounded by spaces: 
        gpr = ' ' + gpr.replace('(', ' ( ').replace(')', ' ) ') + ' '

        
        # translate the spontaneus gene
        if f' {refspont} ' in gpr:
            gpr = gpr.replace(f' {refspont} ', f' spontaneous ')

        
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
    if refspont != 'spontaneous': to_remove = to_remove + [refmodel_t.genes.get_by_id(refspont)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="need to pass in a list")
        cobra.manipulation.remove_genes(refmodel_t, to_remove, remove_reactions=True)
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
            results_presence.append(os.path.exists(f'working/brh/{accession}_align_ref_vs_acc.csv'))
            results_presence.append(os.path.exists(f'working/brh/{accession}_align_acc_vs_ref.csv'))
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



def convert_reference(logger, refmodel, ref_proteome, gene_recovery, refspont):
    
    
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
    if not gene_recovery: edits_dict = None  # load the 'edits_dict'
    else: edits_dict = pickle.load(open('working/rec_broken/edits_dict.pickle', 'rb'))
    create_refgid_to_clusters(logger, refmodel_basename, ref_proteome_basename, edits_dict)
    
    
    # get a copy of the refmodel, and translate its genes to clusters notation. 
    response = translate_refmodel(logger, refmodel, ref_proteome, refspont)
    if response == 1: return 1


    # make traces to keep track of the accessions processed:
    shutil.copyfile('working/annotation/proc_acc.pickle', 'working/brh/proc_acc.pickle')
    
    
    return 0
    