import pickle
import os
import subprocess


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from ..commons import check_cached



def create_combined(logger):
    
    
    # read which proteomes passed the filters: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    proteomes = []  # get the proteomes to parse
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            proteomes.append(proteome)
        
           
    # create foundamental objects:
    sequences_df = []  # future proteins dataframe
    acc_to_seqs = {}  # accession-to-proteins dictionary (1-to-many)
    seq_to_acc = {}  # protein-to-accession dictionary (1-to-1)
             
    
    # parse each proteome:
    sr_list = []
    for proteome in proteomes:
        basename = os.path.basename(proteome)
        accession, _ = os.path.splitext(basename)


        # open the proteome and read its CDSs
        cnt_prots = 0
        with open(proteome, 'r') as cds_handler:
            for seqrecord in SeqIO.parse(cds_handler, "fasta"):
                seqid = seqrecord.id
                seq = seqrecord.seq
                sr = SeqRecord.SeqRecord(seq, id=seqid, description=accession)
                sr_list.append(sr)
                sequences_df.append({'accession': accession, 'cds': seqid, 'aaseq': str(seq)})
                cnt_prots += 1


                # populate the accession-to-proteins dict:
                if accession not in acc_to_seqs.keys(): 
                    acc_to_seqs[accession] = set()
                acc_to_seqs[accession].add(seqid)


                # populate the protein-to-accession dict:
                if seqid not in seq_to_acc.keys(): 
                    seq_to_acc[seqid] = accession
                else: 
                    logger.error(f"The 'seqid' should be unique ('seq_to_acc[{seqid}] = {accession}').")
                    return 1


            # some log messages:
            logger.debug(f"{accession}: {cnt_prots} initial coding sequences found.")


    # write a single fasta with all proteins from all strains
    os.makedirs('working/clustering/', exist_ok=True)
    with open('working/clustering/combined.faa', 'w') as w_handler: 
        count = SeqIO.write(sr_list, w_handler, "fasta")

        
        # write the protein dataframe to file: 
        sequences_df = pnd.DataFrame.from_records(sequences_df)
        sequences_df = sequences_df.set_index('cds', verify_integrity=True)
        sequences_df.to_csv('working/clustering/sequences.csv')
        
        
        # write dictionaries to file:
        with open('working/clustering/acc_to_seqs.pickle', 'wb') as file:
            pickle.dump(acc_to_seqs, file)
        with open('working/clustering/seq_to_acc.pickle', 'wb') as file:
            pickle.dump(seq_to_acc, file)
        
        
        # some log messages:
        logger.debug(f"Total initial coding sequences found: {len(sequences_df)}.")
        
        
        return 0

    

def perform_clustering(logger, cores): 
    
    
    # some log messages: 
    logger.debug("Clustering with " + str(cores) + " cores...")


    # launch the command: 
    with open(f'working/logs/stdout_clustering.txt', 'w') as stdout, open(f'working/logs/stderr_clustering.txt', 'w') as stderr:
        # -d 0: keep the entire sequence ID in the .clstr file
        command = f"""cd-hit \
            -T {cores} -M 0 \
            -i working/clustering/combined.faa \
            -o working/clustering/representatives.faa \
            -g 1 -aL 0.70 -aS 0.70 -c 0.90 -d 0"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
        
    
    # some log messages: 
    logger.debug("Clustering completed.")
    logger.debug("Parsing the clusters...")
    
    
    # create foundamental objects:
    cluster_to_seqs = {}  # cluster to sequences (1-to-many)
    cluster_to_rep = {}  # cluster to representative (1-to-1)
    seq_to_cluster = {}  # sequence to cluster (1-to-1)
    
    
    # parse the clstr_file
    with open('working/clustering/representatives.faa.clstr', 'r') as handler:
        rows = handler.readlines()
        for row in rows: 
            
            
            # the row has a cluster name:
            if row.startswith('>'): 
                cluster = row[1:-1] # '-1' remove the '\n'
                cluster = cluster.replace(' ', '_') # 'Cluster 548' -> 'Cluster_548'
                if cluster not in cluster_to_seqs.keys(): 
                    cluster_to_seqs[cluster] = set()
                    
            
            # the row has a sequence: 
            else: 
                seq, tail = row.split('>',1)[1].rsplit('... ', 1)
                if tail == '*\n':  # if representative:
                    cluster_to_rep[cluster] = seq
                cluster_to_seqs[cluster].add(seq)
                seq_to_cluster[seq] = cluster
                

    # Write the dictionaries to file:
    with open('working/clustering/cluster_to_seqs.pickle', 'wb') as handler:
        pickle.dump(cluster_to_seqs, handler)
    with open('working/clustering/cluster_to_rep.pickle', 'wb') as handler:
        pickle.dump(cluster_to_rep, handler)
    with open('working/clustering/seq_to_cluster.pickle', 'wb') as handler:
        pickle.dump(seq_to_cluster, handler)
        
        
    # load previously created seq_to_acc dictionary: 
    with open('working/clustering/seq_to_acc.pickle', 'rb') as handler:
        seq_to_acc = pickle.load(handler)
        
        
    # create the dictionary rep-to-aaseq:
    rep_to_aaseq = {}
    with open('working/clustering/representatives.faa', 'r') as r_handler: 
        for seqrecord in SeqIO.parse(r_handler, "fasta"):
            seq = seqrecord.seq # <class 'Bio.Seq.Seq'>
            seqid = seqrecord.id
            rep_to_aaseq[seqid] = str(seq)
    with open('working/clustering/rep_to_aaseq.pickle', 'wb') as handler:
        pickle.dump(rep_to_aaseq, handler)
        
        
    # rename the representative sequences:
    with open('working/clustering/representatives.ren.faa', 'w') as w_handler:
        sr_list = []
        with open('working/clustering/representatives.faa', 'r') as r_handler: 
            for seqrecord in SeqIO.parse(r_handler, "fasta"):
                seq = seqrecord.seq # <class 'Bio.Seq.Seq'>
                seqid = seqrecord.id
                new_id = seq_to_cluster[seqid]
                accession = seq_to_acc[seqid]
                sr = SeqRecord.SeqRecord(seq, id=new_id, description=f'{seqid} {accession}')
                sr_list.append(sr)
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
        
    # some log messages: 
    logger.debug("Parsing completed.")
    
    
    return 0



def create_pam(logger):
    
    
    # load previously created dictionaries: 
    with open('working/clustering/cluster_to_seqs.pickle', 'rb') as handler:
        cluster_to_seqs = pickle.load(handler)
    with open('working/clustering/seq_to_acc.pickle', 'rb') as handler:
        seq_to_acc = pickle.load(handler)
    
    
    # create the matrix
    rows = []
    for cluster in cluster_to_seqs.keys():
        
        
        # for each accession, determine the CDSs belonging to this Cluster.
        acc_to_seqs_found = {}
        for seq in cluster_to_seqs[cluster]: 
            accession = seq_to_acc[seq]
            if accession not in acc_to_seqs_found.keys():
                acc_to_seqs_found[accession] = set()
            acc_to_seqs_found[accession].add(seq)
            
        
        # append a row to the dataframe: 
        row = {'cluster': cluster}
        for accession in acc_to_seqs_found.keys():
            row[accession] = ';'.join(acc_to_seqs_found[accession])
        rows.append(row)
        
        
    # write the pam to disk
    pam = pnd.DataFrame.from_records(rows)
    pam = pam.set_index('cluster', verify_integrity=True)
    pam.to_csv('working/clustering/pam.csv')
    
    
    # some log messages: 
    logger.debug("Created an initial presence/absence matrix (PAM) with shape: " + str(pam.shape))
    
    
    # create accession to suffixes dictionary: 
    logger.debug("Storing the suffix for each genome...")
    acc_to_suffix = {}
    for acc in pam.columns:
        suffixes = [i for i in pam[acc].to_list() if type(i) != float]
        suffix = list(set([i.split('_', 1)[0] for i in suffixes]))[0]
        acc_to_suffix[acc] = suffix
    with open('working/clustering/acc_to_suffix.pickle', 'wb') as handler:
        pickle.dump(acc_to_suffix, handler)
        
        
    # create the cluster to relative frequency dictionary: 
    # with the following binary expression, eventual '_stop' are included.
    cluster_to_absfreq = pam.apply(lambda col: col.map(lambda x: 1 if (type(x) != float and x != '') else 0)).sum(axis=1)
    cluster_to_relfreq = round(cluster_to_absfreq / len(pam.columns) * 100, 1)
    cluster_to_relfreq = cluster_to_relfreq.to_dict()
    with open('working/clustering/cluster_to_relfreq.pickle', 'wb') as handler:
        pickle.dump(cluster_to_relfreq, handler)
    
    
    return 0
        


def compute_clusters(logger, cores): 
    
    
    # some log messages:
    logger.info("Clustering the proteins...")
    
    
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/clustering/pam.csv', 
        imp_files=[
            'working/clustering/acc_to_seqs.pickle',
            'working/clustering/cluster_to_rep.pickle',
            'working/clustering/cluster_to_seqs.pickle',
            'working/clustering/seq_to_acc.pickle',
            'working/clustering/seq_to_cluster.pickle',
            'working/clustering/acc_to_suffix.pickle',
            'working/clustering/cluster_to_relfreq.pickle',
            'working/clustering/rep_to_aaseq.pickle',
            'working/clustering/representatives.ren.faa',
            'working/clustering/sequences.csv'])
    if response == 0: return 0
    
    
    # combine together all the sequences from all the strains: 
    response = create_combined(logger)
    if response == 1: return 1
    
    
    # perform the clustering 
    response = perform_clustering(logger, cores)
    if response == 1: return 1
    
    
    # creating the presence/absence matrix (PAM)
    response = create_pam(logger)
    if response == 1: return 1
    
    
    return 0