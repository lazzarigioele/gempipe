import random
import os
import pickle


import pandas as pnd



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



def check_cached(logger, pam_path, imp_files, summary_path=None):
    
    
    # read which proteomes passed the filters: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
        
    # get the accessions to search for: 
    accessions = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            accessions.append(accession)
    accessions = set(accessions)
    
    
    # search for the PAM: 
    if os.path.exists(pam_path):
        pam = pnd.read_csv(pam_path, index_col=0)
        columns = set(list(pam.columns))
        
        
        # search for the optional summary:
        if summary_path != None:
            summary = pnd.read_csv(summary_path, index_col=0)
            rows = set(list(summary.index))
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



def update_pam(logger, pam, module_dir):
    
    
    # define important objects:
    summary = []
    cnt_newgenes = 0
    pam_update = pam.copy()
    
    
    # parse each results file: 
    for file in glob.glob(f'{module_dir}/results/*.csv'):
        accession = file.rsplit('/', 1)[1].replace('.csv', '')
        with open(file, 'r') as r_handler: 
            if r_handler.read() == '""\n':  # if the result csv for this accession is empty: 
                refound_summary.append({'accession': accession, 'n_refound': 0, 'n_stop': 0})
                continue
                
                
        # populate the summary: 
        result = pnd.read_csv(file, sep=',', index_col=0)
        summary.append({
            'accession': accession, 
            'n_refound': len(result[result['gid'].str.contains('_refound')]), 
            'n_stop': len(result[result['gid'].str.contains('_stop')]),
        })

        
        # update the PAM: 
        for cluster in set(result['cluster'].to_list()): 
            new_genes = result[result['cluster']==cluster]['gid'].to_list()
            cnt_newgenes += len(new_genes)
            cell = ';'.join(new_genes)
            pam_update.loc[cluster, accession] = cell
            
        
    # write the summary and update PAM to disk
    summary = pnd.DataFrame.from_records(summary)
    summary = summary.set_index('accession', drop=True, verify_integrity=True)
    summary.to_csv(f'{module_dir}/summary.csv')
    pam_update.to_csv(f'{module_dir}/pam.csv')
    
    
    # write some log messages:
    logger.debug(f'Added {cnt_newgenes} new sequences.')
    
    
    return 0