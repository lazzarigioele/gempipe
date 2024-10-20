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



DEFAULT_TEST_DICT = {
    # amino acids
    'EX_ala__L_e': 'alanine',  # 1
    'EX_arg__L_e': 'arginine',  # 2
    'EX_asp__L_e': 'aspartate',  # 3
    'EX_cys__L_e': 'cysteine' ,  # 4
    'EX_glu__L_e': 'glutamate' ,  # 5
    'EX_gly_e': 'glycine' ,  # 6
    'EX_his__L_e': 'histidine' ,  # 7
    'EX_ile__L_e': 'isoleucine',  # 8
    'EX_leu__L_e': 'leucine',  # 9
    'EX_lys__L_e': 'lysine',  # 10
    'EX_met__L_e': 'methionine',  # 11
    'EX_phe__L_e': 'phenylalanyne',  # 12
    'EX_pro__L_e': 'proline',  # 13
    'EX_ser__L_e': 'serine',  # 14
    'EX_thr__L_e': 'threonine',  # 15
    'EX_trp__L_e': 'tryptophane',  # 16
    'EX_tyr__L_e': 'tyrosine',  # 17
    'EX_val__L_e': 'valine',  # 18
    'EX_asn__L_e': 'asparagine',  # 19
    'EX_gln__L_e': 'glutamine',  # 20
    # vitamins
    'EX_btn_e': 'biotine',  # 1, vitamin B7
    'EX_fol_e': 'folate', # 2, vitamin B9
    'EX_lipoate_e': 'lipoate', # 3, 6,8-Thioctic acid / alpha-Lipoic acid
    'EX_pnto__R_e': 'panthotenate', # 4, vitamin B5
    'EX_pydxn_e': 'pyridoxine',  # 5, form of vitamin B6
    'EX_pydam_e': 'pyridoxamine',  # 6, form of vitamin B6
    'EX_pydx_e': 'pyridoxal',   # form of vitamin B6
    'EX_ribflv_e': 'riboflavin', # 7, vitamin B2
    'EX_thm_e': 'thiamine',  # 8, vitamin B1
    'EX_nac_e': 'nicotinate',  # 9, vitamin PP, vitamin B3, niacin
    'EX_4abz_e': '4_Aminobenzoate', # 10, pABA, vitamin B10
    'EX_cbl1_e': 'cob(I)alamin',   # cobolamine, vitamin B12
    'EX_ascb__L_e': 'ascorbate', # ascorbic acid / vitamin C
}



def auxotropy_simulation(model, seed=False, mode='binary', test_dict=None, model_id=None):
    """
    Function to test auxotrophies in a GSMM. 
    A growth-enabling medium is assumed to be already set up. 
    All compounds -1 in 'test_dict' (aminoacids) will be supplemented.
    
    seed: switch to ModelSEED naming system (not yet impemented)
    mode: 'binary' (1: auxotroph, 0: autotroph) or 'growth': quantitative results from FBA. 
    test_dict:  Dictionary of compounds to test. For example {'EX_ala__L_e': 'alanine', 'EX_arg__L_e': 'arginine', ...}
    model_id: name of the putput column (if None, 'output' will be used)
    """

    
    # get the dictionary of compounds to be tested
    if test_dict == None:
        test_dict = DEFAULT_TEST_DICT
        
    # get the modeled rids: 
    modeled_rids = set([r.id for r in model.reactions])
    
    
    df = [] # list of dict to be converted in pnd dataframe
    if model_id == None: model_id = 'output'
    with model:  # reversible changes. 

        # iterate the compound dictionaries 2 times: 
        # (aa and aa2 are EX_change reactions)
        for aa in test_dict.keys():
            aux_key = f'[aux]{aa[3:-2]}'  # format the dataframe index. For example, from 'EX_glu__L_e' to [aux]glu__L.
            if aa not in modeled_rids:
                df.append({'exchange': aux_key, model_id: None})
                continue 
                
            for aa2 in test_dict.keys():
                if aa2 not in modeled_rids:
                    continue
                    
                if aa2 == aa: 
                    model.reactions.get_by_id(aa2).lower_bound = 0
                else:  # set all other compounds to an arbitrarly high concentration
                    model.reactions.get_by_id(aa2).lower_bound = -1000  # mmol / L

            # perform flux balance analysis. Growth is assumed to be already set as objective. 
            res = model.optimize()

            if res.status == 'optimal' and res.objective_value > 0.001:  # FIVE decimals
                auxotroph = 0 
            else:
                auxotroph = 1

            # save results in a future pnd DataFrame:
            if mode=='binary':
                df.append({'exchange': aux_key, model_id: auxotroph})
            elif mode=='growth':
                if res.status=='optimal': 
                    df.append({'exchange': aux_key, model_id: res.objective_value})
                else: 
                    df.append({'exchange': aux_key, model_id: res.status})
    
    df = pnd.DataFrame.from_records(df)
    df = df.set_index('exchange', drop=True, verify_integrity=True)
    return df




def task_auxotrophy(accession, args):
    
    
    # retrive the arguments:
    outdir = args['outdir']
    skipgf = args['skipgf']
    
    
    # read json/sbml file:
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user asked to skip the strain-specific gapfilling step
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    # perform the simulations: 
    df_results = auxotropy_simulation(ss_model, model_id=accession)
    df_results = df_results.T.reset_index(drop=False).rename(columns={'index': 'accession'})
        
    # it has just 1 row:
    return [df_results.iloc[0].to_dict()]



def strain_auxotrophies_tests(logger, outdir, cores, pam, skipgf):
    
    
    # log some messages
    logger.info("Testing strain-specific auxotrophies for aminoacids and vitamins...")

    
    # check if it's everything pre-computed:
    # not needed!
    
    
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
            itertools.repeat(['accession'] + [i[3:-2] for i in DEFAULT_TEST_DICT.keys()]), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_auxotrophy),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'outdir': outdir, 'skipgf': skipgf}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save the auxotrophyie table in the sae format of 'rpam':
    aux_pam = all_df_combined
    aux_pam = aux_pam.T
    aux_pam.to_csv(outdir + 'aux.csv')