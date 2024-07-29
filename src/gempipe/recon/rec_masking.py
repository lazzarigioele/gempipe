import os
import pickle
import multiprocessing
import itertools
import shutil
import subprocess
import glob


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import check_cached
from ..commons import create_summary
from ..commons import update_pam
from ..commons import extract_aa_seq_from_genome
from ..commons import get_blast_header



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
    with open('working/rec_broken/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords_update = pickle.load(handler)
    logger.debug(f'rec_masking: seq_to_coords: starting from {len(seq_to_coords_update.values())} sequences.')

        
    # now add the new seqs (recovered by this module): 
    for accession in accessions: 
        results_df = pnd.read_csv(f'working/rec_masking/results/{accession}.csv', index_col=0)
        for index, row in results_df.iterrows():
            seq_to_coords_update[row['ID']] = {'accession': row['accession'], 'contig': row['contig'], 'strand': row['strand'], 'start': row['start'], 'end': row['end']}
    logger.debug(f'rec_masking: seq_to_coords: {len(seq_to_coords_update.values())} sequences after the addition of new IDs.')

    
    # save the update dictionary: 
    with open('working/rec_masking/seq_to_coords.pickle', 'wb') as file:
        pickle.dump(seq_to_coords_update, file)



def genome_masking(genome, seq_to_coords):

    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)

    
    # parse each genome:
    sr_list = []
    with open(genome, 'r') as r_handler:                  
        for seqrecord in SeqIO.parse(r_handler, "fasta"):
            contig = seqrecord.id
            seq = seqrecord.seq
            seq_masked = seqrecord.seq  # to be updated
            
            
            # mask the contig from its genes: 
            for seq in seq_to_coords.keys(): 
                if seq_to_coords[seq]['accession'] == accession:
                    if seq_to_coords[seq]['contig'] == contig: 
                        start, end = seq_to_coords[seq]['start'], seq_to_coords[seq]['end']
                        cds_len = end - start  # in this case, 'end' is always greater then 'start'.
                        gene_masked = Seq.Seq(''.join(['N' for i in range(cds_len)]))
                        seq_masked = seq_masked[:start] + gene_masked + seq_masked[end:]
                    
            
            # save the masked squences: 
            sr = SeqRecord.SeqRecord(seq_masked, id=contig, description='')
            sr_list.append(sr)
    with open(f'working/rec_masking/masked_assemblies/{accession}.masked.fna', 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")



def task_recmasking(genome, args):
    
    
    # retrive the arguments:
    pam = args['pam']
    cluster_to_rep = args['cluster_to_rep']
    rep_to_aaseq = args['rep_to_aaseq']
    acc_to_suffix = args['acc_to_suffix']
    seq_to_coords = args['seq_to_coords']
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)
    
    
    # create a query file for each genome: 
    sr_list = []
    with open(f'working/rec_masking/queries/{accession}.query.faa', 'w') as w_handler: 
        for cluster in pam.index:
            cell = pam.loc[cluster, accession]
            if type(cell) == float:  # include only empty clusters
                rep = cluster_to_rep[cluster]
                seq = Seq.Seq(rep_to_aaseq[rep])
                sr = SeqRecord.SeqRecord(seq, id=cluster, description='')
                sr_list.append(sr) 
        count = SeqIO.write(sr_list, w_handler, "fasta")


    # mask the genome from its genes:
    response = genome_masking(genome, seq_to_coords)


    # create a blast database for the genome:
    os.makedirs(f'working/rec_masking/databases/{accession}/', exist_ok=True)
    shutil.copyfile(f'working/rec_masking/masked_assemblies/{accession}.masked.fna', f'working/rec_masking/databases/{accession}/{accession}.masked.fna')  # just the content, not the permissions.
    command = f"""makeblastdb -in working/rec_masking/databases/{accession}/{accession}.masked.fna -dbtype nucl -parse_seqids"""  # '-parse_seqids' is required for 'blastdbcmd'.
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # perform the blast search:
    command = f'''tblastn \
        -query working/rec_masking/queries/{accession}.query.faa \
        -db working/rec_masking/databases/{accession}/{accession}.masked.fna \
        -out working/rec_masking/alignments/{accession}.tsv \
        -outfmt "6 {get_blast_header()}"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # read the alignment: 
    colnames = f'{get_blast_header()}'.split(' ')
    alignment = pnd.read_csv(f'working/rec_masking/alignments/{accession}.tsv', sep='\t', names=colnames )
    alignment['qcov'] = round((alignment['qend'] -  alignment['qstart'] +1)/ alignment['qlen'] * 100, 1)
    alignment['scov'] = round((alignment['send'] -  alignment['sstart'] +1)/ alignment['slen'] * 100, 1)
    
    
    # instantiate key objects: 
    cnt_good_genes = 0
    new_rows = []  # new rows for the sequences_df
    df_result = []  # using the common formatting
    alignment_filtered = []  # future dataframe (blast tbl + seq ID)
    
    
    # group the hsps by cluster:
    cluster_to_indexes = alignment.groupby(['qseqid']).groups
    for cluster in cluster_to_indexes.keys():
        # isolate results for this cluster: 
        indexes = cluster_to_indexes[cluster]
        alignment_cluster = alignment.iloc[indexes, ]
        
        
        # filter using the same thresholds of cd-hit
        alignment_cluster = alignment_cluster[(alignment_cluster['pident'] >= 90) & (alignment_cluster['qcov'] >= 70)]
        
        
        # parse each good hsp: 
        for index, row in alignment_cluster.iterrows():
            cnt_good_genes += 1 
            refound_gid = f"{acc_to_suffix[accession]}_refound_{cnt_good_genes}"


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
                    f'working/rec_masking/databases/{accession}/{accession}.masked.fna',
                    contig, strand, start, end) 
            
            
            # if premature stop, sign the ID
            if len(curr_seq_translated_tostop) / len(curr_seq_translated) < 0.95 :
                refound_gid = refound_gid + '_stop'


            # populate the results dataframe:
            df_result.append({
                'ID': refound_gid, 'cluster': cluster, 
                'accession': accession, 'contig': contig, 'strand': strand, 
                'start': start, 'end': end,
            })
            
            
            # populate the alignment_filtered dataframe
            align_filt_row = row.to_dict()
            align_filt_row['ID'] = refound_gid
            alignment_filtered.append(align_filt_row)
            
            
            # populate the sequences dataframe
            new_rows.append({'cds': refound_gid, 'accession': accession, 'aaseq': str(curr_seq_translated)})


    # save the filtered hsps for this genome:
    alignment_filtered = pnd.DataFrame.from_records(alignment_filtered)
    alignment_filtered.to_csv(f'working/rec_masking/alignments_filtered/{accession}.csv')
    
            
    # save results for this genome:
    df_result = pnd.DataFrame.from_records(df_result)
    df_result.to_csv(f'working/rec_masking/results/{accession}.csv')
    
    
    # return new rows for the sequences_df
    return new_rows



def recovery_masking(logger, cores):
    
    
    # some log messages:
    logger.info("Recovering the genes using genome masking...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/rec_masking/', exist_ok=True)
    os.makedirs('working/rec_masking/queries/', exist_ok=True)
    os.makedirs('working/rec_masking/masked_assemblies/', exist_ok=True)
    os.makedirs('working/rec_masking/databases/', exist_ok=True)
    os.makedirs(f'working/rec_masking/alignments/', exist_ok=True)
    os.makedirs(f'working/rec_masking/alignments_filtered/', exist_ok=True)
    os.makedirs('working/rec_masking/results/', exist_ok=True)
    
    
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_masking/pam.csv',
        summary_path='working/rec_masking/summary.csv',
        imp_files = [
            'working/rec_masking/sequences.csv',
            'working/rec_masking/seq_to_coords.pickle',])
    if response == 0: 
        return 0
    
    
    # load the assets to form the args dictionary:
    pam = pnd.read_csv('working/rec_broken/pam.csv', index_col=0)
    with open('working/clustering/cluster_to_rep.pickle', 'rb') as handler:
        cluster_to_rep = pickle.load(handler)
    with open('working/clustering/rep_to_aaseq.pickle', 'rb') as handler:
        rep_to_aaseq = pickle.load(handler)
    with open('working/clustering/acc_to_suffix.pickle', 'rb') as handler:
        acc_to_suffix = pickle.load(handler)
    with open('working/rec_broken/seq_to_coords.pickle', 'rb') as handler:
        seq_to_coords = pickle.load(handler)


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
            itertools.repeat(task_recmasking),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'pam': pam, 'cluster_to_rep': cluster_to_rep, 'rep_to_aaseq': rep_to_aaseq, 'acc_to_suffix': acc_to_suffix, 'seq_to_coords': seq_to_coords}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # save tabular results:
    sequences_df = pnd.read_csv('working/rec_broken/sequences.csv', index_col=0)
    sequences_df_updated = pnd.concat([sequences_df, all_df_combined], axis=0)
    sequences_df_updated.to_csv('working/rec_masking/sequences.csv')
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # update the pam:
    update_pam(logger, module_dir='working/rec_masking', pam=pam)
    
    
    # update the seq to coordinates dictionary
    update_seq_to_coords(logger)
    
    
    # create the summary:
    create_summary(logger, module_dir='working/rec_masking')
    
    
    return 0
    
    