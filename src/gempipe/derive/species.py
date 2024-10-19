import multiprocessing
import itertools
import os
import shutil


import pandas as pnd
import cobra


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import fba_no_warnings




def task_deriverpam(accession, args):
    
    
    # get the arguments
    outdir = args['outdir']
    skipgf = args['skipgf']
    
    
    # load the model
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user requested to skip the strain-specific gap-filling
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    
    # extract the reactions from the model
    row_dict = {'accession': accession}
    for r in ss_model.reactions:
        row_dict[r.id] = 1

    
    # return the new dataframe row
    return [row_dict]



def derive_rpam(logger, outdir, cores, panmodel, skipgf):

    
    # log some messages
    logger.info("Creating reaction presence/absence matrix...")
    

    # create items for parallelization:
    strains_table = pnd.read_csv(outdir + 'derive_strains.csv', index_col=0)
    items = strains_table.index.to_list()  # get all the accessions
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession'] + [r.id for r in panmodel.reactions]), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_deriverpam),  
            itertools.repeat({'outdir': outdir, 'skipgf': skipgf}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save the rpam
    rpam = all_df_combined
    # now replace missing values with 0. 
    # In 'older' pandas versions, fillna() was doing silent dowcasting (see https://medium.com/@felipecaballero/deciphering-the-cryptic-futurewarning-for-fillna-in-pandas-2-01deb4e411a1)
    # In 'future' pandas versions, silent downcasting will be forbidden. 
    # In 'middle' pandas versions, there is an option to simulate the future behaviour.
    # What follows should guarantee the old behaviour in 'older' and 'middle' pandas versions:
    try:
        with pnd.option_context('future.no_silent_downcasting', True):
            rpam = rpam.fillna(0).infer_objects()   # infer_objects() will do downcasting like in 'older' pandas versions
    except: # OptionError: "No such keys(s): 'future.no_silent_downcasting'" ==> 'older' pandas version
        rpam = rpam.fillna(0)  
    rpam = rpam.astype(int)  # force from float to int.
    rpam = rpam.T  # transpose: reactions as rows.
    rpam.to_csv(outdir + 'rpam.csv')
    
    
    return 0



def task_derivespecies(spp, args):
    
    
    # get the arguments
    strains_table = args['strains_table']
    panmodel = args['panmodel']
    outdir = args['outdir']
    rpam = args['rpam']
    sbml = args['sbml']
    
    
    # filter accessions for this spp
    species_table = strains_table[strains_table['species']==spp]
    
    
    # filter rpam for this spp
    species_rpam = rpam[species_table.index.to_list()]
    
    
    # get rids not present in all the strains of this species:
    to_delete = species_rpam.index[species_rpam.sum(axis=1) != len(species_table)].to_list()
    
    
    # convert rid to panmodel's reaction
    spp_model = panmodel.copy()
    to_delete = [spp_model.reactions.get_by_id(rid) for rid in to_delete]
    
    
    # delete marked reactions
    spp_model.remove_reactions(to_delete)
    
    
    # get some metrics: 
    n_G = len(spp_model.genes)
    n_R = len(spp_model.reactions)
    n_M = len(spp_model.metabolites)
    
    
    # try the FBA: 
    res, obj_value, status = fba_no_warnings(spp_model)
    
    
    # save species specific model to disk
    spp_filename = spp.replace(' ', '_')  # remove spaces for the filename
    cobra.io.save_json_model(spp_model, f'{outdir}/{spp_filename}.json')
    if sbml: cobra.io.write_sbml_model(spp_model, f'{outdir}/{spp_filename}.xml')

    
    # compose the new row
    return [{'species': spp, 'G': n_G, 'R': n_R, 'M': n_M, 'obj_value': obj_value, 'status': status }]



def derive_species_specific(logger, outdir, cores, panmodel, sbml):


    # log some messages
    logger.info("Deriving species-specific models...")
    
   
    # create output dir
    if os.path.exists(outdir + 'species_models/'):
        # always overwriting if already existing
        shutil.rmtree(outdir + 'species_models/')  
    os.makedirs(outdir + 'species_models/', exist_ok=True)
    

    # load the reaction presence/absence table:
    rpam = pnd.read_csv(outdir + 'rpam.csv', index_col=0)
    
    
    # create items for parallelization:
    strains_table = pnd.read_csv(outdir + 'derive_strains.csv', index_col=0)
    items = list(set(strains_table['species'].to_list()))  # get species list
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['species', 'G', 'R', 'M', 'obj_value', 'status']), 
            itertools.repeat('species'), 
            itertools.repeat(logger), 
            itertools.repeat(task_derivespecies),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'strains_table': strains_table, 'panmodel': panmodel, 'outdir': outdir + 'species_models/', 'rpam': rpam, 'sbml': sbml}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save tabular output:
    all_df_combined.to_csv(outdir + 'derive_species.csv')
    
    
    return 0



