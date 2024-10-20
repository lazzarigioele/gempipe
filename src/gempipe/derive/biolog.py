import multiprocessing
import itertools
from importlib import resources
import os
import warnings


import cobra
import pandas as pnd


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import read_refmodel




def biolog_simulation(model, biolog_mappings, seed=False, starting_C='EX_glc__D_e', starting_N='EX_nh4_e', starting_P='EX_pi_e', starting_S='EX_so4_e' ):
    
    
    # get the modeled rids: 
    modeled_rids = set([r.id for r in model.reactions])
    
    
    # iterate each Biolog substrate
    df_results = []  # list of dicts
    for substrate, row in biolog_mappings.iterrows():
        # skip the substrate if 'source'==set():
        for source in eval(row['source']): 
            
            
            # get the 'substituting' and 'substitute' rids:
            exr_after = row['SEED_exchange'] if seed else row['BiGG_exchange']
            if source == 'C': exr_before = starting_C  # 'EX_glc__D_e'
            elif source == 'N': exr_before = starting_N  # 'EX_nh4_e'
            elif source == 'P': exr_before = starting_P  # 'EX_pi_e'
            elif source == 'S': exr_before = starting_S  # 'EX_so4_e'
            else: exr_before = '?'

            # flags for presence of reaction IDs:
            exr_after_present = exr_after in modeled_rids
            exr_before_present = exr_before in modeled_rids
            
            # other columns to fill
            formula = None
            growth = None
            value = None
            status = None
            

            # changes are temporary: 
            with model:
                
                # perform the simulation:
                if exr_before_present:  # the exchange reaction could not be present
                    model.reactions.get_by_id(exr_before).bounds = (0, 1000)
                    
                    # particularly important in case of rich media
                    res_before = model.optimize()
                    objv_before = res_before.objective_value
                    if res_before.status == 'infeasible': objv_before = 0
                    
                    if exr_after_present: # the exchange reaction could not be present
                        r = model.reactions.get_by_id(exr_after)
                        r.bounds = (-1000, 1000)
                        formula = list(r.metabolites)[0].formula
                        
                        # optimize catching the warnings
                        with warnings.catch_warnings(record=True) as w:
                            res = model.optimize()
                        for warning in w:
                            pass
                        
                        value = res.objective_value
                        status = res.status
                        
                        growth = value >= (objv_before + 0.001) and status == 'optimal'
            
            
            # even if they could be simulated, substrate out of the Biolog system
            # are not taken into account in this analysis:
            if eval(row['PM'])==set(): 
                growth = value = status = None
            

            # populate results dictionary:
            df_results.append({
                'substrate': substrate, 'source': source, 'formula': formula,  
                'exr_before': exr_before, 'exr_before_present': exr_before_present,
                'exr_after': exr_after, 'exr_after_present': exr_after_present,
                'growth': growth, 'value': value, 'status': status})
    df_results  = pnd.DataFrame.from_records(df_results)
    
    
    return df_results

    

def task_biolog(accession, args):
    
    
    # retrive the arguments:
    biolog_mappings = args['biolog_mappings']
    outdir = args['outdir']
    skipgf = args['skipgf']
    
    # read json/sbml file:
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user asked to skip the strain-specific gapfilling step
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    # perform the simulations: 
    outfile = outdir + 'biolog/' + accession + '.csv'
    df_results = biolog_simulation(ss_model, biolog_mappings)
    df_results.to_csv(outfile)
    
 
    return [{'accession': accession, 'completed': True}]



def strain_biolog_tests(logger, outdir, cores, pam, panmodel, skipgf):
    
    
    # log some messages
    logger.info("Simulating Biolog's substrates utilization...")
    
    
    # create needed subdirectories:
    os.makedirs(outdir + 'biolog/', exist_ok=True)
    
    
    # check if it's everything pre-computed:
    # not needed!
    
    
    # load the assets to form the args dictionary:
    biolog_mappings = None
    with resources.path('gempipe.assets', 'biolog_mappings.csv' ) as asset_path: 
        biolog_mappings = pnd.read_csv(asset_path, index_col=0)
        
        
    # first of all, test the capabilities of the panmodel: 
    outfile = outdir + 'biolog_panmodel.csv'
    df_results = biolog_simulation(panmodel, biolog_mappings)
    df_results.to_csv(outfile)

    
    # create items for parallelization: 
    items = []
    for accession in pam.columns:
        items.append(accession)
       
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
                          
                          
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession', 'completed']), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_biolog),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'biolog_mappings': biolog_mappings, 'outdir': outdir, 'skipgf': skipgf}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    