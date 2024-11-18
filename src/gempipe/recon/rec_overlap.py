import os
import pickle
import multiprocessing
import itertools
import shutil
import subprocess


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from ..commons import check_cached
from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_blast_header
from ..commons import extract_aa_seq_from_genome
from ..commons import create_summary
from ..commons import update_pam



def update_seq_to_coords(logger): 
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
    
    # get the accessions (passing the quality filters):
    accessions = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            accessions.append(accession)
    

    # create an updateed seq_to_coords dict: 
    with open('working/rec_masking/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords_update = pickle.load(handler)
    logger.debug(f'rec_overlap: seq_to_coords: starting from {len(seq_to_coords_update.values())} sequences.')

        
    # now add the new seqs (recovered by this module): 
    for accession in accessions: 
        results_df = pnd.read_csv(f'working/rec_overlap/results/{accession}.csv', index_col=0)
        for index, row in results_df.iterrows():
            seq_to_coords_update[row['ID']] = {'accession': row['accession'], 'contig': row['contig'], 'strand': row['strand'], 'start': row['start'], 'end': row['end']}
    logger.debug(f'rec_overlap: seq_to_coords: {len(seq_to_coords_update.values())} sequences after the addition of new IDs.')

    
    # save the update dictionary: 
    with open('working/rec_overlap/seq_to_coords.pickle', 'wb') as file:
        pickle.dump(seq_to_coords_update, file)



def task_recoverlap(genome, args):
    
    
    # retrive the arguments:
    pam = args['pam']
    rep_to_aaseq = args['rep_to_aaseq']
    acc_to_suffix = args['acc_to_suffix']
    cluster_to_rep = args['cluster_to_rep']
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)
    
    
    # create a query file for each genome: 
    sr_list = []
    with open(f'working/rec_overlap/queries/{accession}.query.faa', 'w') as w_handler: 
        for cluster in pam.index:
            cell = pam.loc[cluster, accession]
            if type(cell) == float:  # include only empty clusters
                rep = cluster_to_rep[cluster]
                seq = Seq.Seq(rep_to_aaseq[rep])
                sr = SeqRecord.SeqRecord(seq, id=cluster, description='')
                sr_list.append(sr) 
        count = SeqIO.write(sr_list, w_handler, "fasta")
    
    
    # create a blast database for the genome:
    os.makedirs(f'working/rec_overlap/databases/{accession}/', exist_ok=True)
    shutil.copyfile(genome, f'working/rec_overlap/databases/{accession}/{accession}.fna')  # just the content, not the permissions.
    command = f"""makeblastdb -in working/rec_overlap/databases/{accession}/{accession}.fna -dbtype nucl -parse_seqids"""  # '-parse_seqids' is required for 'blastdbcmd'.
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # perform the blast search:
    command = f'''tblastn \
        -query working/rec_overlap/queries/{accession}.query.faa \
        -db working/rec_overlap/databases/{accession}/{accession}.fna \
        -out working/rec_overlap/alignments/{accession}.tsv \
        -outfmt "6 {get_blast_header()}"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # read the masked alignment: 
    colnames = f'{get_blast_header()}'.split(' ')
    masked = pnd.read_csv(f'working/rec_masking/alignments/{accession}.tsv', sep='\t', names=colnames )
    masked['qcov'] = round((masked['qend'] -  masked['qstart'] +1)/ masked['qlen'] * 100, 1)
    masked['scov'] = round((masked['send'] -  masked['sstart'] +1)/ masked['slen'] * 100, 1)
    
    
    # read the NON-masked alignment: 
    colnames = f'{get_blast_header()}'.split(' ')
    not_masked = pnd.read_csv(f'working/rec_overlap/alignments/{accession}.tsv', sep='\t', names=colnames )
    not_masked['qcov'] = round((not_masked['qend'] -  not_masked['qstart'] +1)/ not_masked['qlen'] * 100, 1)
    not_masked['scov'] = round((not_masked['send'] -  not_masked['sstart'] +1)/ not_masked['slen'] * 100, 1)
    

    # focus on genes that are still missing: 
    still_missing = set(not_masked['qseqid'].to_list())
    masked = masked[masked['qseqid'].isin(still_missing)]


    # masked: filter using the same thresholds of cd-hit:
    masked = masked[(masked['pident'] >= 90) & (masked['qcov'] < 70)]
    #masked = masked[masked['evalue'] <= 1e-5]  # if commented, more genes will be recovered
    masked = masked.sort_values('evalue', ascending=True)  # sort by evalue
    masked = masked.reset_index(drop=True)
    
    
    # NON-masked: filter using the same thresholds of cd-hit:
    not_masked = not_masked[(not_masked['pident'] >= 90) & (not_masked['qcov'] >= 70)]
    not_masked = not_masked[not_masked['evalue'] <= 1e-5]
    not_masked = not_masked.sort_values('evalue', ascending=True)
    not_masked = not_masked.reset_index(drop=True)


    # search for sequences that were partially recovered on rec_masked(), 
    # and that now they can be extended to completion. 
    improvements_df = []  # future dataframe showing the elongation to completion. 
    
    
    # group hsps by cluster: 
    groups_not_masked = not_masked.groupby('qseqid').groups
    groups_masked = masked.groupby('qseqid').groups
    for cluster in groups_not_masked.keys():
        if cluster not in groups_masked.keys():
            continue  # impossible to extend
        hpss_not_masked = not_masked.iloc[ groups_not_masked[cluster], ]
        hpss_masked = masked.iloc[ groups_masked[cluster], ]


        # keep the best hsp for this cluster (according to the evalue)
        hpss_masked = hpss_masked.drop_duplicates(subset='qseqid', keep='first')
        hpss_masked = hpss_masked.reset_index(drop=True)
        
        
        # get the contig where to look for a sequence elongation, plus the coverage to beat:
        contig_to_inspect = hpss_masked.loc[0, 'sseqid']
        cov_to_beat = hpss_masked.loc[0, 'qcovhsp']
        
        
        # search for candidates on that contig, having an higher converage:
        hpss_not_masked = hpss_not_masked[hpss_not_masked['sseqid'] == contig_to_inspect]
        hpss_not_masked = hpss_not_masked[hpss_not_masked['qcovhsp'] > cov_to_beat]
        if len(hpss_not_masked)==0: 
            continue  # no candidates found.
        # focus on the best candidated (according to the evalue):
        hpss_not_masked = hpss_not_masked.drop_duplicates(subset='qseqid', keep='first')  

        
        # save the elongation as a two-rows table: 
        elongation_df = pnd.concat([hpss_not_masked, hpss_masked], axis=0)
        elongation_df = elongation_df.reset_index(drop=True)
        elongation_df.loc[1, 'sseqid'] = elongation_df.loc[1, 'sseqid'] + ' (masked)'
        
        
        # check that "start" or "end" of the masked fragment falls between "start" and "end" of the non-masked:
        if elongation_df.loc[1, 'qstart'] < elongation_df.loc[0, 'qstart'] and elongation_df.loc[1, 'qend'] < elongation_df.loc[0, 'qstart']:
            continue
        if elongation_df.loc[1, 'qstart'] > elongation_df.loc[0, 'qend'] and elongation_df.loc[1, 'qend'] > elongation_df.loc[0, 'qend']:
            continue

        # check qlen (length of the representative sequence of the cluster to recover) is at least 100 aa
        if elongation_df.loc[0, 'qlen'] < 100:
            continue

        # check that the fragment (masked contig) is significative (at least 20 aa)
        if elongation_df.loc[1, 'length'] < 20:
            continue

        
        # populate the final improvements dataframe
        improvements_df.append(elongation_df)


    # get and save the final improvements dataframe
    if improvements_df != []:  # if something was recovered
        improvements_df = pnd.concat(improvements_df, axis=0)
        improvements_df = improvements_df.reset_index(drop=True)
    else:  # create an empty dataframe with compatible columns:
        # not_masked and masked should have the same header
        improvements_df = pnd.DataFrame(columns=not_masked.columns)
    improvements_df.to_csv(f'working/rec_overlap/elongations/{accession}.csv')


    # instantiate key objects: 
    cnt_good_genes = 0
    new_rows = []  # new rows for the sequences_df
    df_result = [] # using the common formatting
    
    
    # parse just the elongated sequences: 
    improvements_df = improvements_df[improvements_df['sseqid'].str.contains(' (masked)', regex=False) == False]
    for index, row in improvements_df.iterrows():
        cnt_good_genes += 1
        overlap_gid = f'{acc_to_suffix[accession]}_overlap_{cnt_good_genes}'


        # retrieve the sequence from the genome:
        start = int(row['sstart'])
        end = int(row['send'])
        contig = row["sseqid"]
        contig = contig.split('|')[1]   #split() because blast but some obscure formatting, eg: gb|NIGV01000003.1| for NIGV01000003.1.
        strand = '+'
        if start > end: # if on the other strand, invert the positions. 
            strand = '-'
            start, end = end, start
        curr_seq_translated, curr_seq_translated_tostop = \
            extract_aa_seq_from_genome(
                f'working/rec_overlap/databases/{accession}/{accession}.fna',
                contig, strand, start, end) 


        # if premature stop, sign the ID
        if len(curr_seq_translated_tostop) / len(curr_seq_translated) < 0.95 :
            overlap_gid = overlap_gid + '_stop'


        # populate the results dataframe:
        df_result.append({
            'ID': overlap_gid, 'cluster': row["qseqid"], 
            'accession': accession, 'contig': contig, 'strand': strand, 
            'start': start, 'end': end,
        })
        
        
        # populate the sequences dataframe
        new_rows.append({'cds': overlap_gid, 'accession': accession, 'aaseq': str(curr_seq_translated)})


    # save results for this genome:
    df_result = pnd.DataFrame.from_records(df_result)
    df_result.to_csv(f'working/rec_overlap/results/{accession}.csv')
    
    
    # return new rows for the sequences_df
    return new_rows



def recovery_overlap(logger, cores):
    
    
    # some log messages:
    logger.info("Recovering the genes checking the overlap...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/rec_overlap/', exist_ok=True)
    os.makedirs('working/rec_overlap/queries/', exist_ok=True)
    os.makedirs('working/rec_overlap/databases/', exist_ok=True)
    os.makedirs('working/rec_overlap/alignments/', exist_ok=True)
    os.makedirs('working/rec_overlap/elongations/', exist_ok=True)
    os.makedirs('working/rec_overlap/results/', exist_ok=True)
    
    
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_overlap/pam.csv',
        summary_path='working/rec_overlap/summary.csv',
        imp_files = [
            'working/rec_overlap/sequences.csv',
            'working/rec_overlap/seq_to_coords.pickle',])
    if response == 0: 
        return 0
    
    
    # load the assets to form the args dictionary:
    pam = pnd.read_csv('working/rec_masking/pam.csv', index_col=0)
    with open('working/clustering/acc_to_suffix.pickle', 'rb') as handler:
        acc_to_suffix = pickle.load(handler)
    with open('working/clustering/cluster_to_rep.pickle', 'rb') as handler:
        cluster_to_rep = pickle.load(handler)
    with open('working/clustering/rep_to_aaseq.pickle', 'rb') as handler:
        rep_to_aaseq = pickle.load(handler)
        
        
    # load the previously created species_to_genome: 
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
            itertools.repeat(['cds', 'accession', 'aaseq']), 
            itertools.repeat('cds'), 
            itertools.repeat(logger), 
            itertools.repeat(task_recoverlap),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'pam': pam, 'rep_to_aaseq': rep_to_aaseq, 'acc_to_suffix': acc_to_suffix, 'cluster_to_rep': cluster_to_rep}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # save tabular results:
    sequences_df = pnd.read_csv('working/rec_masking/sequences.csv', index_col=0)
    sequences_df_updated = pnd.concat([sequences_df, all_df_combined], axis=0)
    sequences_df_updated.to_csv('working/rec_overlap/sequences.csv')
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # update the pam:
    update_pam(logger, module_dir='working/rec_overlap', pam=pam)
    
    
    # update the seq to coordinates dictionary
    update_seq_to_coords(logger)
    
    
    # create the summary:
    create_summary(logger, module_dir='working/rec_overlap')
    
    
    return 0