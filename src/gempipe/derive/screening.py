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
    
    return 0
    
    
    
def get_sources_by_class(model):
    
    sources_by_class = {'C': set(), 'N': set(), 'P': set(), 'S': set()}
    for r in model.reactions: 
        if len(r.metabolites)==1 and list(r.metabolites)[0].id.endswith('_e'):
            m = list(r.metabolites)[0]
            formula = m.formula
            # avoid confusion with 'C':
            formula = formula.replace('Ca', '').replace('Co', '').replace('Cu', '').replace('Cd', '').replace('Cr', '').replace('Cs', '').replace('Cl', '')   
            # avoid confusion with 'N':
            formula = formula.replace('Na', '').replace('Nb', '').replace('Ni', '').replace('Ne', '')
            # avoid confusion with 'P':
            formula = formula.replace('Pd', '').replace('Pt', '').replace('Pb', '').replace('Po', '')
            # avoid confusion with 'S':
            formula = formula.replace('Sc', '').replace('Si', '').replace('Sn', '').replace('Sb', '').replace('Se', '')
            
            if 'C' in formula: sources_by_class['C'].add(r.id)
            if 'N' in formula: sources_by_class['N'].add(r.id)
            if 'P' in formula: sources_by_class['P'].add(r.id)
            if 'S' in formula: sources_by_class['S'].add(r.id)
    
    return sources_by_class



def cnps_simulation(model, seed=False, mode='binary', sources_by_class=None, model_id=None, starting_C='EX_glc__D_e', starting_N='EX_nh4_e', starting_P='EX_pi_e', starting_S='EX_so4_e'):
    """
    Function to test utilization of C-N-P-S substrates in a GSMM. 
    A growth-enabling medium is assumed to be already set up. 
    
    seed: switch to ModelSEED naming system (not yet impemented)
    mode: 'binary' (1: auxotroph, 0: autotroph) or 'growth': quantitative results from FBA. 
    sources_by_class:  Dictionary of compounds to test. For example {'C': {'EX_ala__L_e', ...}, 'N': {'EX_ala__L_e', ...}}
    model_id: name of the putput column (if None, 'output' will be used)
    """    
    
    # get the dictionary of compounds to be tested
    if sources_by_class == None:
        sources_by_class = get_sources_by_class(model)
        
    # get the modeled rids: 
    modeled_rids = set([r.id for r in model.reactions])
    
    
    df = [] # list of dict to be converted in pnd dataframe
    if model_id == None: model_id = 'output'
    
    
    for sub_class, starting in zip(['C','N','P','S'], [starting_C, starting_N, starting_P, starting_S]):
        if starting in modeled_rids:   # For example, 
            for exr_after in sources_by_class[sub_class]:
                
                with model:  # reversible changes 
                    # close the original substrate
                    model.reactions.get_by_id(starting).lower_bound = 0
                    
                    # first FBA to be later compared:
                    res_before = model.optimize()
                
                    # open the alternative substrate:
                    model.reactions.get_by_id(exr_after).lower_bound = -1000
                    
                    # second FBA for camparison:
                    res_after = model.optimize()
                    
                    if res_before.status=='optimal' and res_after.status=='optimal' and res_after.objective_value >= (res_before.objective_value + 0.001):
                        can_use = 1
                    else:
                        can_use = 0
                        
                    # save results in a future pnd DataFrame:
                    sub_key = f'[{sub_class}]{exr_after[3:-2]}'
                    if mode=='binary':
                        df.append({'exchange': sub_key, model_id: can_use})
                    elif mode=='growth':
                        if res_after.status=='optimal': 
                            df.append({'exchange': sub_key, model_id: res_after.objective_value})
                        else: 
                            df.append({'exchange': sub_key, model_id: res_after.status})

    
    df = pnd.DataFrame.from_records(df)
    df = df.set_index('exchange', drop=True, verify_integrity=True)
    return df
    
    
    
def task_cnps(accession, args):
    
    
    # retrive the arguments:
    outdir = args['outdir']
    skipgf = args['skipgf']
    sources_by_class = args['sources_by_class']
    
    
    # read json/sbml file:
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user asked to skip the strain-specific gapfilling step
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    # perform the simulations: 
    df_results = cnps_simulation(ss_model, model_id=accession, sources_by_class=sources_by_class)
    df_results = df_results.T.reset_index(drop=False).rename(columns={'index': 'accession'})
        
    # it has just 1 row:
    return [df_results.iloc[0].to_dict()]

    
    
def strain_cnps_tests(logger, outdir, cores, pam, panmodel, skipgf):
    
    
    # log some messages
    logger.info("Testing strain-specific consumption of C-N-P-S substrates...")
    sources_by_class = get_sources_by_class(panmodel)   
    
    
    # ge the header for the results table:
    header = []
    for sub_class in sources_by_class.keys():
        for sub in sources_by_class[sub_class]:
            if sub_class == 'C': sub_key = f'[C]{sub[3:-2]}'
            if sub_class == 'N': sub_key = f'[N]{sub[3:-2]}'
            if sub_class == 'P': sub_key = f'[P]{sub[3:-2]}'
            if sub_class == 'S': sub_key = f'[S]{sub[3:-2]}'
            header.append(sub_key)

    
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
            itertools.repeat(['accession'] + header), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_cnps),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'outdir': outdir, 'skipgf': skipgf, 'sources_by_class': sources_by_class}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save the auxotrophyie table in the sae format of 'rpam':
    cnps_pam = all_df_combined
    cnps_pam = cnps_pam.T
    cnps_pam.to_csv(outdir + 'cnps.csv')
    
    return 0