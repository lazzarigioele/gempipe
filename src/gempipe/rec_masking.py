import os
import pickle
import multiprocessing
import itertools
import shutil
import subprocess
import glob


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from .commons import chunkize_items
from .commons import load_the_worker
from .commons import gather_results
from .commons import check_cached
from .commons import update_pam
from .commons import extract_aa_seq_from_genome




def genome_masking(genome):
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)


    ### TODO
        
    
    # parse each genome:
    sr_list = []
    with open(genome, 'r') as r_handler:                  
        for seqrecord in SeqIO.parse(r_handler, "fasta"):
            contig = seqrecord.id
            seq = seqrecord.seq
            seq_masked = seqrecord.seq  # to be updated
            
            
            # mask the contig from its genes: 
            if contig in contig_to_gff.keys():  # a contig may have no CDS.
                for gene_coords in contig_to_gff[contig]: 
                    start, end = gene_coords[0], gene_coords[1]
                    cds_len = end - start  # in this case, 'end' is always greater then 'start'. 
                    gene_masked = Seq.Seq(''.join(['N' for i in range(cds_len)]))
                    seq_masked = seq_masked[:start] + gene_masked + seq_masked[end:]
                    
            
            # save the masked squences: 
            sr = SeqRecord.SeqRecord(seq_masked, id=contig, description='')
            sr_list.append(sr)
    with open(f'working/rec_masking/masked_assemblies/{accession}.masked.fna', 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
    
    return 0



def task_recmasking(genome, args):
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)
    
    
    # create a query file for each genome: 
    sr_list = []
    with open(f'working/rec_masking/queries/{accession}.query.faa', 'w') as w_handler: 
        for cluster in args['pam'].index:
            cell = args['pam'].loc[cluster, accession]
            if type(cell) == float: 
                rep = args['cluster_to_rep'][cluster]
                seq = Seq.Seq(args['sequences_df'].loc[rep, 'aaseq'])
                sr = SeqRecord.SeqRecord(seq, id=cluster, description='')
                sr_list.append(sr) 
        count = SeqIO.write(sr_list, w_handler, "fasta")


    # mask the genome from its genes:
    response = genome_masking(genome)


    # create a blast database for the genome:
    os.makedirs(f'working/rec_masking/databases/{accession}/', exist_ok=True)
    # shutil.copy is specified to copy permission bits, while shutil.copyfile is just to copy the file contents.
    shutil.copyfile(f'working/rec_masking/masked_assemblies/{accession}.masked.fna', f'working/rec_masking/databases/{accession}/{accession}.masked.fna')  
    # '-parse_seqids' is required for later use of 'blastdbcmd'.
    command = f"""makeblastdb -in working/rec_masking/databases/{accession}/{accession}.masked.fna -dbtype nucl -parse_seqids"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # perform the blast search:
    command = f'''tblastn \
        -query working/rec_masking/queries/{accession}.query.faa \
        -db working/rec_masking/databases/{accession}/{accession}.masked.fna \
        -out working/rec_masking/alignments/{accession}.tsv \
        -outfmt "6 qseqid sseqid pident ppos length qlen slen qstart qend sstart send evalue bitscore qcovhsp"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    # select the best hits: 
    cnt_good_genes = 0
    new_rows = []
    df_result = []
    colnames = 'qseqid sseqid pident ppos length qlen slen qstart qend sstart send evalue bitscore qcovhsp'.split(' ')
    alignment = pnd.read_csv(f'working/rec_masking/alignments/{accession}.tsv', sep='\t', names=colnames )
    cluster_to_indexes = alignment.groupby(['qseqid']).groups
    for cluster in cluster_to_indexes.keys():
        # isolate results for this cluster: 
        indexes = cluster_to_indexes[cluster]
        alignment_cluster = alignment.iloc[indexes, ]
        
        
        # filter using the same thresholds of cd-hit
        alignment_cluster = alignment_cluster[(alignment_cluster['qcovhsp'] >= 70) & (alignment_cluster['pident'] >= 90)]
        for index, row in alignment_cluster.iterrows():
            cnt_good_genes += 1 
            refound_gid = f"{args['acc_to_suffix'][accession]}_refound_{cnt_good_genes}"


            # retrieve the sequence from the genome:
            start = int(row['sstart'])
            end = int(row['send'])
            contig = row["sseqid"]
            #contig = contig.split('|')[1]   #split() because blast but some obscure formatting, eg: gb|NIGV01000003.1| for NIGV01000003.1.
            strand = '+'
            if start > end: # if on the other strand, invert the positions. 
                strand = '-'
                start, end = end, start
            curr_seq_translated, curr_seq_translated_tostop = \
                extract_aa_seq_from_genome(
                    f'working/rec_masking/databases/{accession}/{accession}.masked.fna',
                    contig, strand, start, end) # TODO
            if len(curr_seq_translated_tostop) / len(curr_seq_translated) < 0.95 :
                refound_gid = refound_gid + '_stop'


            # finally save all results:
            df_result.append({
                'gid': refound_gid, 'cluster': cluster, 
                'accession': accession, 'contig': row['sseqid'], 
                'start': row['sstart'], 'end': row['send'],
                'coverage': row['qcovhsp'], 'pident': row['pident'], 'ppos': row['ppos'],
                'evalue': row['evalue'], 'bitscore': row['bitscore'], 
            })
            new_rows.append({'cds': refound_gid, 'accession': accession, 'aaseq': str(curr_seq_translated_tostop)})


    # save results for this genome:
    df_result = pnd.DataFrame.from_records(df_result)
    df_result.to_csv(f'working/rec_masking/results/{accession}.csv')
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
    os.makedirs('working/rec_masking/results/', exist_ok=True)
    
    
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_masking/pam.csv',
        summary_path='working/rec_masking/summary.csv',
        imp_files = ['working/rec_masking/sequences.csv'])
    if response == 0: return 0
    
    
    # load the assets ot form the args dictionary:
    pam = pnd.read_csv('working/rec_broken/pam.csv', index_col=0)
    sequences_df = pnd.read_csv('working/rec_broken/sequences.csv', index_col=0)
    with open('working/rec_broken/cluster_to_rep.pickle', 'rb') as handler:
        cluster_to_rep = pickle.load(handler)
    with open('working/rec_broken/acc_to_suffix.pickle', 'rb') as handler:
        acc_to_suffix = pickle.load(handler)


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
            itertools.repeat(task_recmasking),
            itertools.repeat({'pam': pam, 'cluster_to_rep': cluster_to_rep, 'sequences_df': sequences_df, 'acc_to_suffix': acc_to_suffix,}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # save tabular results:
    sequences_df = pnd.concat([sequences_df, all_df_combined], axis=0)
    sequences_df.to_csv('working/rec_masking/sequences.csv')
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join()  
    
    
    # update the pam and get the summary:
    update_pam(logger, module_dir='working/rec_masking', pam=pam)
    
    
    return 0
    
    