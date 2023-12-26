import os
import pickle
import subprocess


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq



def func_annot(logger, cores, outdir): 
    
    
    # create subdirs without overwriting
    os.makedirs('working/annotation/', exist_ok=True)
    os.makedirs('working/annotation/db/', exist_ok=True)
    
    
    # print some log messages
    logger.info('Performing functional annotation of the representative sequences...')
    
    
    # load the final PAM: 
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    
    
    # check if all the output where already computed: 
    clusters = list(pam.index)
    if os.path.exists('working/annotation/pan.emapper.annotations'):
        if os.path.exists('working/annotation/representatives.faa'):
            seq_ids = [] 
            with open('working/annotation/representatives.faa', 'r') as r_handler:
                for seqrecord in SeqIO.parse(r_handler, "fasta"):
                    seq_ids.append(seqrecord.id)
            if set(seq_ids) == set(clusters):
                # log some message: 
                logger.info('Found all the needed files already computed. Skipping this step.')
                # signal to skip this module:
                return 0
    
    
    # load the sequences resources: 
    with open('working/clustering/cluster_to_rep.pickle', 'rb') as handler:
        cluster_to_rep = pickle.load(handler)
    with open('working/clustering/rep_to_aaseq.pickle', 'rb') as handler:
        rep_to_aaseq = pickle.load(handler)
        
        
    # parse the pam to create a single input fasta files with representative sequences: 
    sr_list = []
    for cluster in pam.index:
        rep = cluster_to_rep[cluster]
        aaseq = Seq.Seq(cluster_to_rep[rep])
        sr = SeqRecord.SeqRecord(aaseq, id=cluster, description=f'({rep})')
        sr_list.append(sr)
    with open(f'working/annotation/representatives.faa', 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
        
    # check if the database already exists:
    if (not os.path.exists('working/annotation/db/eggnog_proteins.dmnd')) or (not os.path.exists('working/annotation/db/eggnog.db')):
        logger.info("The database for functional annotation is missing. It will be dowloaded now...")
        with open(f'working/logs/stdout_funcdownload.txt', 'w') as stdout, open(f'working/logs/stderr_funcdownload.txt', 'w') as stderr: 
            command = f"""download_eggnog_data.py -y --data_dir working/annotation/db/"""
            process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
            process.wait()
        logger.info("Download completed. Now executing the annotation...")
        
    
    # execute the command.
    # --dbmem loads the whole eggnog.db sqlite3 annotation database during the annotation step, and 
    # therefore requires ~44 GB of memory. It is recommanded when annotating a large number of sequences.
    # download_eggnog_data.py : This will download the eggNOG annotation database (along with the taxa databases), 
    # and the database of eggNOG proteins for Diamond searches.
    with open(f'working/logs/stdout_funcannot.txt', 'w') as stdout, open(f'working/logs/stderr_funcannot.txt', 'w') as stderr: 
        command = f"""emapper.py \
            --cpu {cores} \
            --override \
            --data_dir working/annotation/db/ \
            -i working/annotation/representatives.faa \
            -m diamond \
            --itype proteins \
            --trans_table 11 \
            --excel \
            --output pan; mv pan.* working/annotation/"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()

    
    
    return 0