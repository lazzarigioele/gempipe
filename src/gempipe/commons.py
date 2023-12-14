import random
import os


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
        new_row = function(item, args)
        df_combined.append(new_row)
        

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