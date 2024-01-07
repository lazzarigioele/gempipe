import random
import os
import pickle
import glob
import subprocess
import hashlib



import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq



def chunkize_items(items, cores):
    
    
    # divide items in chunks: 
    random.shuffle(items)  # randomly re-order items
    nitems_inchunk = int(len(items) / cores)
    if len(items) % cores !=0: nitems_inchunk += 1
    chunks = [items[x *nitems_inchunk : (x+1)* nitems_inchunk] for x in range(cores)]
    
    
    return chunks



def load_the_worker(arguments):
        
        
    # get the arguments:
    items = arguments[0]
    worker = arguments[1] +1   
    columns = arguments[2]
    index = arguments[3]
    logger = arguments[4]
    function = arguments[5]
    args = arguments[6]


    # iterate over each item: 
    df_combined = []
    cnt_items_processed = 0
    for item in items:


        # perform the annotation for this genome: 
        new_rows = function(item, args)
        df_combined = df_combined + new_rows
        

        # notify the logging process: 
        cnt_items_processed += 1
        logger.debug(f"W#{worker}-PID {os.getpid()}: {round(cnt_items_processed/len(items)*100, 1)}%")


    # join the tabular results of each item:
    if df_combined == []:  # this worker was started empty.
        df_combined = pnd.DataFrame(columns = columns)
    else: df_combined = pnd.DataFrame.from_records(df_combined)
    df_combined = df_combined.set_index(index, verify_integrity=True)
    return df_combined



def gather_results(results):
    
    
    # perform final concatenation of the tabular results:
    all_df_combined = []
    for result in results: 
        if isinstance(result, pnd.DataFrame):
            all_df_combined.append(result)
    all_df_combined = pnd.concat(all_df_combined, axis=0)
    
    
    return all_df_combined



def get_retained_accessions():
    # to be called after the genomes filtering.
    
    
    accessions = set()
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        for species in species_to_proteome.keys(): 
            for proteome in species_to_proteome[species]:
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)
                accessions.add(accession)
    return accessions



def check_cached(logger, pam_path, imp_files, summary_path=None):
    
    
    # get the accessions retained:
    accessions = get_retained_accessions()
    
    
    # search for the PAM: 
    if os.path.exists(pam_path):
        pam = pnd.read_csv(pam_path, index_col=0)
        columns = set(list(pam.columns))
        
        
        # search for the optional summary:
        if summary_path != None:
            if os.path.exists(summary_path): 
                summary = pnd.read_csv(summary_path, index_col=0)
                rows = set(list(summary.index))
            else: return None
        else: rows = columns
            
            
        # check if accessions are the same (no less, no more):
        if accessions == columns == rows:
            
            
            # check the presence of important files:
            if all([os.path.exists(i) for i in imp_files]):
                # log some message: 
                logger.info('Found all the needed files already computed. Skipping this step.')
                
                
                # signal to skip this module:
                return 0
                
            
    return None



def create_summary(logger, module_dir):
    
    
    # get the accessions retained:
    accessions = get_retained_accessions()
    
    
    # parse each results file: 
    summary = []
    for file in glob.glob(f'{module_dir}/results/*.csv'):
        accession = file.rsplit('/', 1)[1].replace('.csv', '')
        if accession not in accessions: 
            continue  # other files present from previous runs.
        with open(file, 'r') as r_handler: 
            if r_handler.read() == '""\n':  # if the result csv for this accession is empty: 
                refound_summary.append({'accession': accession, 'n_refound': 0, 'n_frag': 0, 'n_overlap': 0, 'n_stop': 0})
                continue
                
                
        # populate the summary: 
        result = pnd.read_csv(file, sep=',', index_col=0)
        summary.append({
            'accession': accession, 
            'n_refound': len(result[result['ID'].str.contains('_refound')]), 
            'n_frag': len(result[result['ID'].str.contains('_frag')]), 
            'n_overlap': len(result[result['ID'].str.contains('_overlap')]), 
            'n_stop': len(result[result['ID'].str.contains('_stop')]),
        })

        
    # write the summary to disk
    summary = pnd.DataFrame.from_records(summary)
    summary = summary.set_index('accession', drop=True, verify_integrity=True)
    summary.to_csv(f'{module_dir}/summary.csv')
    
    
    return 0



def update_pam(logger, module_dir, pam):
    
    
    # get the accessions retained:
    accessions = get_retained_accessions()
    
    
    # define important objects:
    cnt_newgenes = 0
    pam_update = pam.copy()
    
    
    # parse each results file: 
    for file in glob.glob(f'{module_dir}/results/*.csv'):
        accession = file.rsplit('/', 1)[1].replace('.csv', '')
        if accession not in accessions: 
            continue  # other files present from previous runs.
        with open(file, 'r') as r_handler: 
            if r_handler.read() == '""\n':  # if the result csv for this accession is empty: 
                continue
        # populate the summary: 
        result = pnd.read_csv(file, sep=',', index_col=0)


        # update the PAM: 
        for cluster in set(result['cluster'].to_list()): 
            new_genes = result[result['cluster']==cluster]['ID'].to_list()
            cnt_newgenes += len(new_genes)
            cell = ';'.join(new_genes)
            pam_update.loc[cluster, accession] = cell
            
        
    # write the updated PAM to disk
    pam_update.to_csv(f'{module_dir}/pam.csv')
    
    
    # write some log messages:
    logger.debug(f'Added {cnt_newgenes} new sequences.')
    
    
    return 0



def extract_aa_seq_from_genome(db, contig, strand, start, end ):
    
    
    # execute the command
    command = f'''blastdbcmd -db {db} -entry "{contig}" -range "{start}-{end}" -outfmt "%s"'''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    process.wait()
    
    
    # read the output from stdout
    curr_stream = process.stdout.read()
    curr_stream = curr_stream.decode('utf8')  # See https://stackoverflow.com/q/41918836
    curr_stream = curr_stream.rstrip()  # remove end of line
    
    
    # reverse complement if on the negative strand
    seq = Seq.Seq(curr_stream)
    if strand == '-': 
        seq = seq.reverse_complement()

        
    # trim the sequences to make it multiple of three, otherwise I get the following warning: 
    # BiopythonWarning: Partial codon, len(sequence) not a multiple of three. 
    seq_trimmed = seq[:len(seq) - (len(seq) % 3)]


    # translate to stop:
    seq_translated = seq_trimmed.translate()
    seq_translated_tostop = seq_trimmed.translate(to_stop=True)
    
    
    return seq_translated, seq_translated_tostop



def get_blast_header():
    
    
    # to standardize all the blast subprocesses
    return "qseqid sseqid pident ppos length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp"


def get_md5_string(filepath):
    
    
    with open(filepath, 'rb') as file:  # 'rb' good also for txt files.
        md5 = hashlib.md5()
        md5.update(file.read())  # warning: no chunks; keep attention with large files.
        md5_string = md5.hexdigest()
    return md5_string