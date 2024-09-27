import os
import shutil 
import subprocess
import pickle
import multiprocessing
import itertools


import pandas as pnd


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import check_cached
from ..commons import create_summary
from ..commons import extract_aa_seq_from_genome
from ..commons import get_blast_header



def alignment_to_couples(accession, cluster_to_relfreq, seq_to_cluster): 
    
    
    # read the alignment with extra columns: 
    colnames = f'{get_blast_header()}'.split(' ')
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


    # group hsps by cluster:
    df_couples = []
    groups = alignment.groupby('sseqid').groups
    for cluster in groups.keys():
        hsps = alignment.iloc[ groups[cluster], ]
        # if a CDS is repeated for the same cluster, take the best hit.
        hsps = hsps.drop_duplicates(subset='qseqid', keep='first')
        # get the length of the representative seq (it's constant)
        slen = hsps['slen'].values[0] 


        # search for couples having progressive numbers:
        prognums = list(hsps['prognum'].values)
        for i in prognums:
            if i in blacklist: continue
            if i+1 in prognums and i+2 not in prognums and i-1 not in prognums:
                # get the two pieces of this protein: 
                good_couple = hsps[hsps['prognum'].isin([i, i+1])]   
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
                relfreq_cluster = cluster_to_relfreq[cluster]
                relfreq_piece1 = cluster_to_relfreq[seq_to_cluster[good_couple['qseqid'].values[0]]]
                relfreq_piece2 = cluster_to_relfreq[seq_to_cluster[good_couple['qseqid'].values[1]]]

                
                # if this couple respect the thresholds: 
                if  overall_cov >= 70 and \
                    sup <= 30 and \
                    rqlen1 < 90 and rqlen2 < 90 and \
                    relfreq_cluster > relfreq_piece1 and relfreq_cluster > relfreq_piece2:
                    
                    
                    # include in results and update the blacklist:
                    df_couples.append(good_couple)
                    blacklist.add(i)
                    blacklist.add(i+1)

                    
                # consider just the first progressive couple (that is, with the highest evalues):
                break  


    # write results to disk: 
    if df_couples == []:  # create an empty dataframe
        df_couples = [pnd.DataFrame(columns=alignment.columns)]
    df_couples = pnd.concat(df_couples, axis=0)
    df_couples = df_couples.drop('prognum', axis=1)
    df_couples = df_couples.reset_index(drop=True)
    df_couples.to_csv(f'working/rec_broken/couples/{accession}.csv')
    
    
    return df_couples


    
def get_updated_column(pam, accession, df_couples, cluster_to_relfreq, seq_to_cluster, seq_to_coords):
    
    
    # this module will release a new updated pam.
    # the new pam will be created gluing together singular columns.
    # columns will be first translated (.T) to be compliant with the .commons lib.
    pam_column = pam.loc[: , [accession]].copy()
    
    
    # save all the new mappings
    edits_dict = {}
    
    
    # parse each couple to trace the jumping of protein pieces:
    with open(f'working/rec_broken/edits/{accession}.txt', "w") as w_handler:
        groups = df_couples.groupby('sseqid').groups
        for cluster in groups.keys():
            couple = df_couples.iloc[ groups[cluster], ]


            # convert the original cell of the gained cluster to set
            ori_cell = pam_column.loc[cluster, accession]
            if type(ori_cell)==float: ori_cell_set = set()  # empty cell
            else: ori_cell_set = set(ori_cell.split(';'))


            # get the two pieces (prefix and progressive number)
            cds_ids = couple['qseqid'].to_list()
            cds_prognums = sorted([i.rsplit('_', 1)[1] for i in cds_ids]) 
            prefix = cds_ids[0].split('_', 1)[0]
            
            
            # determine if the two CDSs are coming from the same contig AND the same strand.
            # WIth the following implementation, broken proteins are recovered only if they come from the same contig. 
            # Further controls are applied during populate_results_df(), raising errors.
            if seq_to_coords[cds_ids[0]]['contig'] != seq_to_coords[cds_ids[1]]['contig']:
                continue
            if seq_to_coords[cds_ids[0]]['strand'] != seq_to_coords[cds_ids[1]]['strand']:
                continue


            # remove eventual pieces from the gained cluster cell, then add the glued protein
            ori_cell_depleted = ori_cell_set - set(cds_ids)
            new_cell = set([f'{prefix}_frag_' + '_'.join(cds_prognums)]).union(ori_cell_depleted)
            new_cell = ';'.join(new_cell)

            
            # get the freq of the gained cluster, log the change, and apply to pam
            rel_freq = cluster_to_relfreq[cluster]
            print(f'{cluster} ({rel_freq}%): {ori_cell} --> {new_cell}', file=w_handler)
            pam_column.loc[cluster, accession] = new_cell


            # now update the cells of the two pieces, following the same logic above.
            for cds in cds_ids: 
                cluster2 = seq_to_cluster[cds]
                ori_cell2 = pam_column.loc[cluster2, accession]
                ori_cell2_set = set() if type(ori_cell2)==float else set(ori_cell2.split(';'))
                new_cell2 = set(ori_cell2_set) - set(cds_ids) 
                new_cell2 = ';'.join(new_cell2)
                rel_freq2 = cluster_to_relfreq[cluster2]
                print(f'\t{cluster2} ({rel_freq2}%): {ori_cell2} --> {new_cell2}', file=w_handler)
                pam_column.loc[cluster2, accession] = new_cell2
                
                
            # new: create a new seq_to_cluster dictionary
            for cds in cds_ids: 
                edits_dict[f'{cds}:{seq_to_cluster[cds]}'] = f'{prefix}_frag_' + '_'.join(cds_prognums) + f':{cluster}'
    
    
    # save the dictionary: 
    with open(f'working/rec_broken/edits/{accession}.pickle', 'wb') as file:
        pickle.dump(edits_dict, file)
    
    
    return pam_column[accession]



def task_recbroken(genome, args):
    
    
    # retrive the arguments:
    pam = args['pam']
    cluster_to_relfreq = args['cluster_to_relfreq']
    seq_to_cluster = args['seq_to_cluster']
    seq_to_coords = args['seq_to_coords']
    # WARNING: seq_to_coords can be really heavy (eg 200 MB) when genomes are hundreds.
    # This can significanlty slow down the creation of child processess, and RAM comsumption can be really high.
    
    
    # get the accession and proteome file:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)
    proteome = f'working/proteomes/{accession}.faa'
    
    
    # create a database for later extraction of recovered sequences: 
    os.makedirs(f'working/rec_broken/databases/{accession}/', exist_ok=True)
    shutil.copyfile(genome, f'working/rec_broken/databases/{accession}/{accession}.fna')  # just the content, not the permissions.
    command = f"""makeblastdb -in working/rec_broken/databases/{accession}/{accession}.fna -dbtype nucl -parse_seqids""" # '-parse_seqids' is required for 'blastdbcmd'.
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # perform the blastp: proteins on representatives:  
    command = f'''blastp \
        -query {proteome} \
        -db working/rec_broken/representatives/representatives.ren.faa \
        -out working/rec_broken/alignments/{accession}.tsv \
        -outfmt "6 {get_blast_header()}"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # parse the alignment to get the good couples.
    df_couples = alignment_to_couples(accession, cluster_to_relfreq, seq_to_cluster)
    
    
    # parse the good couples to update the pam column:
    pam_column = get_updated_column(pam, accession, df_couples, cluster_to_relfreq, seq_to_cluster, seq_to_coords)

    
    # return new rows for load_the_worker():
    row = pam_column.to_dict()
    row['accession'] = accession
    return [row]
    
    

def populate_results_df(logger):
    
    
    # load the previously created dictionaries: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    with open('working/coordinates/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords = pickle.load(handler)
        
    
    # parse each accession:
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            results_df = []
            
            
            # get the couple
            df_couples = pnd.read_csv(f'working/rec_broken/couples/{accession}.csv', index_col=0)
            groups = df_couples.groupby('sseqid').groups
            for cluster in groups.keys():
                couple = df_couples.iloc[ groups[cluster], ]
                
                
                # get the new sequence id
                cds_ids = couple['qseqid'].to_list()
                cds_prognums = sorted([i.rsplit('_', 1)[1] for i in cds_ids]) 
                prefix = cds_ids[0].split('_', 1)[0]
                new_seq_id = f'{prefix}_frag_' + '_'.join(cds_prognums)
                
                
                # get the accession
                couple_accessions = set([seq_to_coords[seq]['accession'] for seq in cds_ids])
                if len(couple_accessions) != 1:
                    logger.error(f"Found different accessions in this couple: {cds_ids}.")
                    return 1
                couple_accession = list(couple_accessions)[0]
                
                # repeating controls for contigs and strands, first time in get_updated_column(),
                # is needed because, while the generation of working/rec_broken/edits can take into account contigs/strands,
                # the generation of working/rec_broken/couples (used here) cannot. 
                
                # get the contig
                couple_contigs = set([seq_to_coords[seq]['contig'] for seq in cds_ids])
                if len(couple_contigs) != 1:  # same control implemented in get_updated_column()
                    #logger.debug(f"Found different contigs in this couple: {cds_ids} (accession {couple_accession}).")
                    #return 1
                    continue
                couple_contig = list(couple_contigs)[0]
                
                
                # get the strand
                couple_strands = set([seq_to_coords[seq]['strand'] for seq in cds_ids])
                if len(couple_strands) != 1:  # same control implemented in get_updated_column()
                    #logger.debug(f"Found different strands in this couple: {cds_ids} (accession {couple_accession}).")
                    #return 1
                    continue
                couple_strand = list(couple_strands)[0]
                
                
                # get the start - end
                start_1 = seq_to_coords[cds_ids[0]]['start']
                start_2 = seq_to_coords[cds_ids[1]]['start']
                end_1 = seq_to_coords[cds_ids[0]]['end']
                end_2 = seq_to_coords[cds_ids[1]]['end']
                if start_1 > end_1: 
                    logger.info(f"Found start_1 < end_1 in this couple: {cds_ids} (accession {couple_accession}).")
                    return 1
                if start_2 > end_2: 
                    logger.info(f"Found start_2 < end_2 in this couple: {cds_ids} (accession {couple_accession}).")
                    return 1
                couple_start = min([start_1, start_2])
                couple_end = max([end_1, end_2])
                
                
                # write the results dataframe
                results_df.append({
                    'ID': new_seq_id, 'cluster': cluster, 
                    'accession': couple_accession, 'contig': couple_contig, 'strand': couple_strand, 
                    'start': couple_start, 'end': couple_end
                })
            results_df = pnd.DataFrame.from_records(results_df)
            results_df.to_csv(f'working/rec_broken/results/{accession}.csv')
    
    
    return 0
    
    

def update_seq_to_coords(logger): 
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
    
    # get the good accessions (passing the quality filters):
    good_accessions = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            good_accessions.append(accession)
    
    
    # parse the couples/log files to get the seqs ID to erease:
    to_erease = []
    for accession in good_accessions: 
        df_couples = pnd.read_csv(f'working/rec_broken/couples/{accession}.csv', index_col=0)
        to_erease = to_erease + df_couples['qseqid'].to_list()
    
        
    # create an updateed seq_to_coords dict: 
    with open('working/coordinates/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords = pickle.load(handler)
    logger.debug(f'rec_broken: seq_to_coords: starting from {len(seq_to_coords.values())} sequences.')
    seq_to_coords_update = {}
        
    
    # remove seqs from the updated dictionary:
    for seq in seq_to_coords.keys(): 
        attribs = seq_to_coords[seq] 
        if attribs['accession'] not in good_accessions: 
            continue  # remove all seq IDs belonging to accessions not passing the quality filter. 
        if seq in to_erease: 
            continue  # remove seq IDs belonging to the good couples discovered by this module. 
        seq_to_coords_update[seq] = attribs
    logger.debug(f'rec_broken: seq_to_coords: {len(seq_to_coords_update.values())} after removing filtered accessions and frag couples.')

        
    # now add the new seqs (recovered by this module): 
    for accession in good_accessions: 
        results_df = pnd.read_csv(f'working/rec_broken/results/{accession}.csv', index_col=0)
        for index, row in results_df.iterrows():
            seq_to_coords_update[row['ID']] = {'accession': row['accession'], 'contig': row['contig'], 'strand': row['strand'], 'start': row['start'], 'end': row['end']}
    logger.debug(f'rec_broken: seq_to_coords: {len(seq_to_coords_update.values())} sequences after the addition of new IDs.')

    
    # save the update dictionary: 
    with open('working/rec_broken/seq_to_coords.pickle', 'wb') as file:
        pickle.dump(seq_to_coords_update, file)
        
        
        
def join_edits_dict(logger):
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
    
    # get the good accessions (passing the quality filters):
    good_accessions = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            good_accessions.append(accession)
    
    
    # parse the couples/log files to get the seqs ID to erease:
    edits_dict = {}
    for accession in good_accessions: 
        with open(f'working/rec_broken/edits/{accession}.pickle', 'rb') as handler:
            curr_edits = pickle.load(handler)
        for key, value in curr_edits.items():
            edits_dict[key] = value
            
            
    # save the dictionary: 
    with open(f'working/rec_broken/edits_dict.pickle', 'wb') as file:
        pickle.dump(edits_dict, file)
        
        

def update_sequences(logger): 
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
    
    # get the accessions to work with: 
    accessions = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            accessions.append(accession)
            
    
    # parse the couples/log files to get the seqs ID to erease:
    to_erease = []
    for accession in accessions: 
        df_couples = pnd.read_csv(f'working/rec_broken/couples/{accession}.csv', index_col=0)
        to_erease = to_erease + df_couples['qseqid'].to_list()
            
            
    # create an updated df, removing the seqs to erease:
    sequences_df = pnd.read_csv('working/clustering/sequences.csv' , index_col=0)
    logger.debug(f'rec_broken: sequences dataframe: starting from {len(sequences_df)} sequences.')
    sequences_df_updated = sequences_df.copy()
    sequences_df_updated = sequences_df_updated.drop(to_erease)
    logger.debug(f'rec_broken: sequences dataframe: {len(sequences_df_updated)} after removing frag couples.')
    
    
    # now add the new seqs
    new_rows = []
    for accession in accessions: 
        results_df = pnd.read_csv(f'working/rec_broken/results/{accession}.csv', index_col=0)
        for index, row in results_df.iterrows():
            contig = row['contig']
            strand = row['strand']
            start = row['start']
            end = row['end']
            seq, seq_tostop = extract_aa_seq_from_genome(
                f'working/rec_broken/databases/{accession}/{accession}.fna', 
                contig, strand, start, end)
            new_rows.append({'cds': row['ID'], 'accession': accession, 'aaseq': seq})
    new_rows = pnd.DataFrame.from_records(new_rows)
    new_rows = new_rows.set_index('cds', drop=True, verify_integrity=True)
    sequences_df_updated = pnd.concat([sequences_df_updated, new_rows])
    logger.debug(f'rec_broken: seq_to_coords: {len(sequences_df_updated)} sequences after the addition of new IDs.')

    
    # save the update version:
    sequences_df_updated.to_csv('working/rec_broken/sequences.csv')
    


def recovery_broken(logger, cores):
    
    
    # some log messages:
    logger.info("Recovering the proteins broken in two pieces...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/rec_broken/', exist_ok=True)
    os.makedirs('working/rec_broken/edits/', exist_ok=True)
    os.makedirs('working/rec_broken/representatives/', exist_ok=True)
    os.makedirs('working/rec_broken/databases/', exist_ok=True)
    os.makedirs(f'working/rec_broken/alignments/', exist_ok=True)
    os.makedirs('working/rec_broken/couples/', exist_ok=True)
    os.makedirs('working/rec_broken/results/', exist_ok=True)
    
    
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_broken/pam.csv',
        summary_path='working/rec_broken/summary.csv',
        imp_files = [
            'working/rec_broken/sequences.csv',
            'working/rec_broken/seq_to_coords.pickle',
            'working/rec_broken/edits_dict.pickle',])
    if response == 0: 
        return 0
    

    # copy representative sequences (all) and make a database
    shutil.copyfile(f'working/clustering/representatives.ren.faa', f'working/rec_broken/representatives/representatives.ren.faa') 
    command = f"""makeblastdb -in working/rec_broken/representatives/representatives.ren.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # load the assets to form the args dictionary:
    pam = pnd.read_csv('working/clustering/pam.csv', index_col=0)
    with open('working/clustering/cluster_to_relfreq.pickle', 'rb') as handler:
        cluster_to_relfreq = pickle.load(handler)
    with open('working/clustering/seq_to_cluster.pickle', 'rb') as handler:
        seq_to_cluster = pickle.load(handler)
    with open('working/coordinates/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords = pickle.load(handler)
    
    
    # load the previously created species_to_proteome: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)
        
        
    # create items for parallelization: 
    items = []
    for species in species_to_genome.keys(): 
        for genome in species_to_genome[species]: 
            items.append(genome)
            
            
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
            itertools.repeat(task_recbroken),  # will return a new updated pam.
            itertools.repeat({'pam': pam, 'cluster_to_relfreq': cluster_to_relfreq, 'seq_to_cluster': seq_to_cluster, 'seq_to_coords': seq_to_coords}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # get the updated pam, and remove empty rows:
    pam_updated = all_df_combined.T
    pam_updated = pam_updated.replace({'': None})
    empty_rows_bool = pam_updated.isnull().all(axis=1)
    empty_rows_df = pam_updated.loc[empty_rows_bool]
    logger.debug(f"rec_broken: ended up with {len(empty_rows_df)} devoided clusters.")
    pam_updated = pam_updated.drop(empty_rows_df.index)
    pam_updated.to_csv('working/rec_broken/pam.csv')
    
    
    # get the results dataframe following conventions
    response = populate_results_df(logger)
    if response == 1: return 1

    
    # update the squence to coordinates dict (removing filtered genomes, and protein frags)
    update_seq_to_coords(logger)  # creates 'working/rec_broken/seq_to_coords.pickle'
    
    
    # update the sequences dataframe (removing protein frags)
    update_sequences(logger)  # creates 'working/rec_broken/sequences.csv'
    
    
    # join together all the strain specific edits dicts
    join_edits_dict(logger)  # creates 'working/rec_broken/edits_dict.pickle'
    
    
    # create a summary for this module reading the results dataframes: 
    create_summary(logger, module_dir='working/rec_broken/')
    

    
    return 0