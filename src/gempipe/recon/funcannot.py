import os
import pickle
import subprocess


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from ..commons import get_outdir



def func_annot(logger, cores, outdir, dbs, dbmem): 
    
    
    # create subdirs without overwriting
    os.makedirs('working/annotation/', exist_ok=True)
    dbs = get_outdir(dbs)  # append a '/' and create if necessary.
    
 
    
    # print some log messages
    logger.info('Performing functional annotation of the representative sequences...')
    
    
    # load the final PAM: 
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    
    
    # get the accessions retained:
    accessions = set()
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        for species in species_to_proteome.keys(): 
            for proteome in species_to_proteome[species]:
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)
                accessions.add(accession)
    
    
    # check if all the output where already computed:
    if os.path.exists('working/annotation/proc_acc.pickle'):
        with open('working/annotation/proc_acc.pickle', 'rb') as handler:
            proc_acc = pickle.load(handler) 
        if accessions == proc_acc == set(list(pam.columns)):
            if os.path.exists('working/annotation/pan.emapper.annotations'):
                if os.path.exists('working/annotation/representatives.faa'):
                    seq_ids = [] 
                    with open('working/annotation/representatives.faa', 'r') as r_handler:
                        for seqrecord in SeqIO.parse(r_handler, "fasta"):
                            seq_ids.append(seqrecord.id)
                    clusters = list(pam.index)
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
        aaseq = Seq.Seq(rep_to_aaseq[rep])
        sr = SeqRecord.SeqRecord(aaseq, id=cluster, description=f'({rep})')
        sr_list.append(sr)
    with open(f'working/annotation/representatives.faa', 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
        
    # check if the database already exists:
    if (not os.path.exists(dbs + 'eggnog_proteins.dmnd')) or (not os.path.exists(dbs + 'eggnog.db')):
        logger.info("The database for functional annotation is missing. It will be dowloaded now...")
        with open(f'working/logs/stdout_funcdownload.txt', 'w') as stdout, open(f'working/logs/stderr_funcdownload.txt', 'w') as stderr: 
            command = f"""download_eggnog_data.py -y --data_dir {dbs}"""
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
            --data_dir {dbs} \
            -i working/annotation/representatives.faa \
            -m diamond \
            --itype proteins \
            --trans_table 11 \
            --excel {'--dbmem' if dbmem else ''} \
            --output pan; mv pan.* working/annotation/"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()

        
    # make traces to keep track of the accessions processed
    proc_acc = set(list(pam.columns))
    with open('working/annotation/proc_acc.pickle', 'wb') as handler:
        pickle.dump(proc_acc, handler)
    
    
    return 0