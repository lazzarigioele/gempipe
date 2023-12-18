import os
import shutil 
import subprocess
import pickle
import multiprocessing
import itertools


import pandas as pnd


from .commons import chunkize_items
from .commons import load_the_worker
from .commons import gather_results



def task_recbroken(proteome, args):
    
    
    # get the basename without extension:
    basename = os.path.basename(proteome)
    accession, _ = os.path.splitext(basename)
    
    
    ###
    ### TMP 
    ### TODO
    ###
    if not os.path.exists(f'working/rec_broken/alignments/{accession}.tsv'):
        # perform the blastp 
        command = f'''blastp \
            -query {proteome} \
            -db working/rec_broken/database/representatives.ren.faa \
            -out working/rec_broken/alignments/{accession}.tsv \
            -outfmt "6 qseqid sseqid pident ppos length qlen slen qstart qend sstart send evalue bitscore qcovhsp"
        '''
        process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        process.wait()
        
    
    # this module will release a new updated pam.
    # the new pam will be created guing together singular columns.
    # columns will be first translated (.T) to be compliant with the .commons lib.
    pam_column = args['pam'].loc[: , [accession]].copy()
    
    
    # read the alignment with extra columns: 
    colnames = 'qseqid sseqid pident ppos length qlen slen qstart qend sstart send evalue bitscore qcovhsp'.split(' ')
    alignment = pnd.read_csv(f'working/rec_broken/alignments/{accession}.tsv', sep='\t', names=colnames )
    alignment['qcov'] = round((alignment['qend'] -  alignment['qstart'] +1)/ alignment['qlen'] * 100, 1)
    alignment['scov'] = round((alignment['send'] -  alignment['sstart'] +1)/ alignment['slen'] * 100, 1)
    
    
    # filter and sort by evalue (sorting is important because we consider each CDS just 1 time)
    alignment = alignment[alignment['pident'] >= 90]
    alignment = alignment[alignment['qcov'] >= 70]
    alignment = alignment[alignment['evalue'] <= 1e-5]
    alignment = alignment.sort_values('evalue', ascending=True)
    alignment = alignment.reset_index(drop=True)


    # create the 'prognum' (progressive number) column:
    alignment['prognum'] = None
    for index, row in alignment.iterrows(): 
        alignment.loc[index, 'prognum'] = int(row['qseqid'].split('_', 1)[1])


    # create a blacklist because each CDS must be taken just 1 time: 
    blacklist = set()


    # detect proteins broken in two pieces:
    df_result = []
    groups = alignment.groupby('sseqid').groups
    for cluster in groups.keys():
        hpss = alignment.iloc[ groups[cluster], ]
        # if a CDS is repeated for the same cluster, take the best hit.
        hpss = hpss.drop_duplicates(subset='qseqid', keep='first')
        # get the length of the representative seq (it's constant)
        slen = hpss['slen'].values[0] 


        # search for couples having progressive numbers:
        prognums = list(hpss['prognum'].values)
        for i in prognums:
            if i in blacklist: continue
            if i+1 in prognums and i+2 not in prognums and i-1 not in prognums:
                # get the two pieces of this protein: 
                good_couple = hpss[hpss['prognum'].isin([i, i+1])]   
                # below we compute several metrics: 

                
                # overall coverage % (include the gap between the two pieces):
                overall_cov = ( max(good_couple['send']) - min(good_couple['sstart']) +1 ) / slen * 100


                # superimposition % between the two pieces: 
                if min(good_couple['send']) >= max(good_couple['sstart']):
                    sup = ( min(good_couple['send']) - max(good_couple['sstart']) +1 ) / slen * 100
                else: sup = 0


                # compute the gap % between the two pieces:
                if min(good_couple['send']) < max(good_couple['sstart']):
                    gap = ( max(good_couple['sstart']) - min(good_couple['send']) -1 ) / slen * 100
                else: gap = 0


                # compute the query len relative to the subject (constant)
                rqlen1 = good_couple['qlen'].values[0] / slen * 100
                rqlen2 = good_couple['qlen'].values[1] / slen * 100


                # get the relative frequencies:
                relfreq_cluster = args['cluster_to_relfreq'][cluster]
                relfreq_piece1 = args['cluster_to_relfreq'][args['seq_to_cluster'][good_couple['qseqid'].values[0]]]
                relfreq_piece2 = args['cluster_to_relfreq'][args['seq_to_cluster'][good_couple['qseqid'].values[1]]]

                
                # if this couple respect the thresholds: 
                if  overall_cov >= 70 and \
                    sup <= 30 and \
                    rqlen1 < 90 and rqlen2 < 90 and \
                    relfreq_cluster > relfreq_piece1 and relfreq_cluster > relfreq_piece2:
                    
                    
                    # include in results and update the blacklist:
                    df_result.append(good_couple)
                    blacklist.add(i)
                    blacklist.add(i+1)

                    
                # consider just the first progressive couple (that is, with the highest evalues):
                break  


    # write results to disk: 
    df_result = pnd.concat(df_result, axis=0)
    df_result = df_result.drop('prognum', axis=1)
    df_result = df_result.reset_index(drop=True)
    df_result.to_csv(f'working/rec_broken/results/{accession}.tsv', sep='\t')

    
    # trace the jumping of protein pieces:
    with open(f'working/rec_broken/editlogs/{accession}.txt', "w") as w_handler:
        groups = df_result.groupby('sseqid').groups
        for cluster in groups.keys():
            couple = df_result.iloc[ groups[cluster], ]


            # PHASE 1: update the cell for this clustes.

            ori_cell = pam_column.loc[cluster, accession]
            if type(ori_cell)==float: ori_cell_set = set()  # empty cell
            else: ori_cell_set = set(ori_cell.split(';'))


            cds_ids = couple['qseqid'].to_list()
            cds_prognums = [i.rsplit('_', 1)[1] for i in cds_ids] 
            prokka_prefix = cds_ids[0].split('_', 1)[0]


            ori_cell_depleted = ori_cell_set - set(cds_ids)
            new_cell = set([f'{prokka_prefix}_frag_' + '_'.join(cds_prognums)]).union(ori_cell_depleted)
            new_cell = ';'.join(new_cell)

            
            rel_freq = args['cluster_to_relfreq'][cluster]
            print(f'{cluster} ({rel_freq}%): {ori_cell} --> {new_cell}', file=w_handler)
            pam_column.loc[cluster, accession] = new_cell


            # PHASE 2: update the cells of the originally involvold clusters.

            for cds in cds_ids: 
                cluster_to_edit = args['seq_to_cluster'][cds]
                ori_cell2 = pam_column.loc[cluster_to_edit, accession]
                ori_cell2_set = set() if type(ori_cell2)==float else set(ori_cell2.split(';'))
                new_cell2 = set(ori_cell2_set) - set(cds_ids) 
                new_cell2 = ';'.join(new_cell2)

                
                rel_freq = args['cluster_to_relfreq'][cluster_to_edit]
                print(f'\t{cluster} ({rel_freq}%): {ori_cell2} --> {new_cell2}', file=w_handler)
                pam_column.loc[cluster_to_edit, accession] = new_cell2

    
    # return new rows for load_the_worker():
    row = pam_column[accession].to_dict()
    row['accession'] = accession
    return [row]
    
    
    

def recovery_broken(logger, cores):
    
    
    # some log messages:
    logger.info("Recovering the proteins broken in two pieces...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/rec_broken/', exist_ok=True)
    os.makedirs('working/rec_broken/editlogs/', exist_ok=True)
    os.makedirs('working/rec_broken/database/', exist_ok=True)
    os.makedirs(f'working/rec_broken/alignments/', exist_ok=True)
    os.makedirs('working/rec_broken/results/', exist_ok=True)
    
    """
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_broken/pam.csv',
        summary_path='working/rec_broken/summary.csv',
        imp_files = ['working/rec_broken/sequences.csv'])
    if response == 0: return 0
    """

    # copy representative sequences (all) and make a database
    shutil.copyfile(f'working/clustering/representatives.ren.faa', f'working/rec_broken/database/representatives.ren.faa') 
    command = f"""makeblastdb -in working/rec_broken/database/representatives.ren.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # load the assets to form the args dictionary:
    pam = pnd.read_csv('working/clustering/pam.csv', index_col=0)
    with open('working/clustering/seq_to_cluster.pickle', 'rb') as handler:
        seq_to_cluster = pickle.load(handler)
    
    
    # get cluster to relative frequencies dict: 
    # with the following binary expression, eventual '_stop' are included.
    cluster_to_absfreq = pam.applymap(lambda x: 1 if (type(x) != float and x != '') else 0 ).sum(axis=1)
    cluster_to_relfreq = round(cluster_to_absfreq / len(pam.columns) * 100, 1)
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
        
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
            itertools.repeat(['accession'] + list(pam.T.columns)), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_recbroken),
            itertools.repeat({'pam': pam, 'cluster_to_relfreq': cluster_to_relfreq, 'seq_to_cluster': seq_to_cluster, }),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    """
    # get the updated pam:
    update_pam(logger, pam, module_dir='working/rec_masking')
    """
    pam_updated = all_df_combined.T
    pam_updated.to_csv('working/rec_broken/pam.csv')
    
    
    return 0