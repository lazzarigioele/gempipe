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
from ..commons import fba_no_warnings


from ..interface.medium import reset_growth_env
from ..interface.gaps import can_synth
from ..interface.gaps import get_objectives



DEFAULT_TEST_DICT = {
    # amino acids
    'ala__L': 'alanine',  # 1
    'arg__L': 'arginine',  # 2
    'asp__L': 'aspartate',  # 3
    'cys__L': 'cysteine' ,  # 4
    'glu__L': 'glutamate' ,  # 5
    'gly': 'glycine' ,  # 6
    'his__L': 'histidine' ,  # 7
    'ile__L': 'isoleucine',  # 8
    'leu__L': 'leucine',  # 9
    'lys__L': 'lysine',  # 10
    'met__L': 'methionine',  # 11
    'phe__L': 'phenylalanyne',  # 12
    'pro__L': 'proline',  # 13
    'ser__L': 'serine',  # 14
    'thr__L': 'threonine',  # 15
    'trp__L': 'tryptophane',  # 16
    'tyr__L': 'tyrosine',  # 17
    'val__L': 'valine',  # 18
    'asn__L': 'asparagine',  # 19
    'gln__L': 'glutamine',  # 20
    # vitamins
    'btn': 'biotine',  # 1, vitamin B7
    #'fol': 'folate', # 2, vitamin B9
    'thf': 'tetrahydrofolate', # active form of B9
    'lipoate': 'lipoate', # 3, 6,8-Thioctic acid / alpha-Lipoic acid
    'pnto__R': 'panthotenate', # 4, vitamin B5
    'pydxn': 'pyridoxine',  # 5, form of vitamin B6
    'pydam': 'pyridoxamine',  # 6, form of vitamin B6
    'pydx': 'pyridoxal',   # form of vitamin B6
    'ribflv': 'riboflavin', # 7, vitamin B2
    'thm': 'thiamine',  # 8, vitamin B1
    'nac': 'nicotinate',  # 9, vitamin PP, vitamin B3, niacin
    '4abz': '4_Aminobenzoate', # 10, pABA, vitamin B10
    'cbl1': 'cob(I)alamin',   # cobolamine, vitamin B12
    'ascb__L': 'ascorbate', # ascorbic acid / vitamin C
}



def auxotropy_simulation(model, seed=False, mode='binary', test_dict=None, model_id=None):
    """
    Function to test auxotrophies in a GSMM. 
    A growth-enabling medium is assumed to be already set up. 
    All compounds -1 in 'test_dict' (aminoacids) will be supplemented.
    
    seed: switch to ModelSEED naming system (not yet impemented)
    mode: 'binary' (1: auxotroph, 0: prototroph) or 'growth': quantitative results from FBA. 
    test_dict:  Dictionary of compounds to test. For example {'EX_ala__L_e': 'alanine', 'EX_arg__L_e': 'arginine', ...}
    model_id: name of the putput column (if None, 'output' will be used)
    """
    
    # implementation in versions <= v1.37.5 was based on the assumption that all the 
    # aminoacids/vitamins had their exchange reaction. From v1.37.6 it is sufficient that 
    # a metabolite appears in the cytosol to test for its production. 

    
    # get the dictionary of compounds to be tested
    if test_dict == None:
        test_dict = DEFAULT_TEST_DICT
        
    # get the modeled rids / mids: 
    modeled_rids = set([r.id for r in model.reactions])
    modeled_mids = set([m.id for m in model.metabolites])
    
    
    df = [] # list of dict to be converted in pnd dataframe
    if model_id == None: model_id = 'output'
    with model:  # reversible changes. 

        # iterate the compound dictionaries 2 times: 
        # (aa and aa2 are EX_change reactions)
        for aa in test_dict.keys():
            aux_key = f'[aux]{aa}'  # format the dataframe index. For example, from 'EX_glu__L_e' to '[aux]glu__L'
            aa_c = f'{aa}_c'   # For example, from 'EX_glu__L_e' to 'glu__L_c'.
            
            if aa_c not in modeled_mids:
                auxotroph = 1
                obj_value = 0.0
                status = 'missing'
            
            else:    
                for aa2 in test_dict.keys():
                    aa2_c = f'{aa2}_c'   # For example, from 'EX_glu__L_e' to 'glu__L_c'.
                    if aa2_c not in modeled_mids:
                        continue

                    EX_aa2 = f'EX_{aa2}_e'
                    if aa2_c == aa_c: 
                        if EX_aa2 in modeled_rids: 
                            model.reactions.get_by_id(EX_aa2).lower_bound = 0
                    else:  # set all other compounds to an arbitrarly high concentration
                        if EX_aa2 in modeled_rids: 
                            model.reactions.get_by_id(EX_aa2).lower_bound = -1000  # mmol / L

                # perform flux balance analysis
                binary, obj_value, status = can_synth(model, aa_c)

                if status == 'optimal' and obj_value > 0.00001:  # FIVE decimals
                    auxotroph = 0 
                else:
                    auxotroph = 1

            # save results in a future pnd DataFrame:
            if mode=='binary':
                df.append({'exchange': aux_key, model_id: auxotroph})
            elif mode=='growth':
                if res.status=='optimal': 
                    df.append({'exchange': aux_key, model_id: obj_value})
                else: 
                    df.append({'exchange': aux_key, model_id: status})
    
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
            itertools.repeat(['accession'] + [i for i in DEFAULT_TEST_DICT.keys()]), 
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
            if formula == None: 
                continue
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



def cnps_simulation(model, seed=False, mode='binary', sources_by_class=None, model_id=None, starting_C='EX_glc__D_e', starting_N='EX_nh4_e', starting_P='EX_pi_e', starting_S='EX_so4_e', cnps_minmed=0):
    """
    Function to test utilization of C-N-P-S substrates in a GSMM. 
    A growth-enabling medium is assumed to be already set up. 
    
    seed: switch to ModelSEED naming system (not yet impemented)
    mode: 'binary' (1: can grow, 0: cannot grow) or 'growth': quantitative results from FBA. 
    sources_by_class:  Dictionary of compounds to test. For example {'C': {'EX_ala__L_e', ...}, 'N': {'EX_ala__L_e', ...}}
    model_id: name of the putput column (if None, 'output' will be used)
    cnps_minmed: base the analysis on a minimal medium yielding at least the specified objective value.
        If ``False`` user-defined medium will be applied. 
    """    

    
    # get the dictionary of compounds to be tested
    if sources_by_class == None:
        sources_by_class = get_sources_by_class(model)
        
        
    # if autostarting, a minimal medium is applied and starting C, N, P and,S sources are automatically defined:
    if cnps_minmed != 0.0: 
        
        
        # create a backup for the medium:
        medium_backup = {}
        for r in model.reactions:
            if len(r.metabolites)==1 and list(r.metabolites)[0].id.endswith('_e'):
                medium_backup[r.id] = r.bounds
                
        
        # define the minumum medium: 
        min_medium = cobra.medium.minimal_medium(model, cnps_minmed, minimize_components=True)
        if min_medium is None:  # usually associated with "Minimization of medium was infeasible."
            min_medium = {}  # future 'pandas.core.series.Series'
            for r in model.reactions:   # convert the current medium to a pnd.Series
                if len(r.metabolites)==1 and list(r.metabolites)[0].id.endswith('_e'):
                    if r.lower_bound < 0:
                        min_medium[r.id] = r.lower_bound * -1.0  # positive sign
            min_medium = pnd.Series(min_medium)
        
        
        # apply the min medium: 
        min_medium = min_medium.sort_values(ascending=False)
        reset_growth_env(model)
        for exr_id, lb in min_medium.items():
            model.reactions.get_by_id(exr_id).lower_bound = -lb
            
            
        # define the starting sources:
        for exr_recommended, sub_class in zip(['EX_glc__D_e', 'EX_nh4_e', 'EX_pi_e', 'EX_so4_e'], ['C', 'N', 'P', 'S']):
            if exr_recommended in list(model.medium.keys()):
                if sub_class == 'C': starting_C = exr_recommended
                if sub_class == 'N': starting_N = exr_recommended
                if sub_class == 'P': starting_P = exr_recommended
                if sub_class == 'S': starting_S = exr_recommended
            else:
                for exr_id in list(model.medium.keys()):
                    if exr_id in sources_by_class[sub_class]:
                        if sub_class == 'C': 
                            starting_C = exr_id
                            break
                        if sub_class == 'N': 
                            starting_N = exr_id
                            break
                        if sub_class == 'P': 
                            starting_P = exr_id
                            break
                        if sub_class == 'S': 
                            starting_S = exr_id
                            break
    
            
        
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
                    res_before, obj_value_before, status_before = fba_no_warnings(model)
                
                    # open the alternative substrate:
                    model.reactions.get_by_id(exr_after).lower_bound = -1000
                    
                    # second FBA for camparison:
                    res_after, obj_value_after, status_after = fba_no_warnings(model)
                    
                    if status_before=='infeasible' and status_after=='optimal' and obj_value_after >= 0.001: 
                        can_use = 1
                    elif status_after=='infeasible':
                        can_use = 0
                    elif status_before=='optimal' and status_after=='optimal' and obj_value_after >= (obj_value_before + 0.001):
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

    
    
    # restore medium from backup 
    if cnps_minmed != 0.0:
        for r in model.reactions:
            if len(r.metabolites)==1 and list(r.metabolites)[0].id.endswith('_e'):
                r.bounds = medium_backup[r.id]
    
    
    
    df = pnd.DataFrame.from_records(df)
    df = df.set_index('exchange', drop=True, verify_integrity=True)
    return df
    
    
    
def task_cnps(accession, args):
    
    
    # retrive the arguments:
    outdir = args['outdir']
    skipgf = args['skipgf']
    sources_by_class = args['sources_by_class']
    cnps_minmed = args['cnps_minmed']
    
    
    # read json/sbml file:
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user asked to skip the strain-specific gapfilling step
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    # perform the simulations: 
    df_results = cnps_simulation(ss_model, model_id=accession, sources_by_class=sources_by_class, cnps_minmed=cnps_minmed)
    df_results = df_results.T.reset_index(drop=False).rename(columns={'index': 'accession'})
        
    # it has just 1 row:
    return [df_results.iloc[0].to_dict()]

    
    
def strain_cnps_tests(logger, outdir, cores, pam, panmodel, skipgf, cnps_minmed):
    
    
    # log some messages
    logger.info("Testing strain-specific consumption of C-N-P-S substrates...")
    if cnps_minmed != 0.0: 
        logger.info(f"A minimal medium leading to objective value >= {cnps_minmed} will be used for each strain.")
    
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
            itertools.repeat({'outdir': outdir, 'skipgf': skipgf, 'sources_by_class': sources_by_class, 'cnps_minmed': cnps_minmed}),
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



def biosynth_simulation(model, model_id=None, biosynth=0.5):
    """
    Function to test biosynthesis of metabolites in a model. 
    A growth-enabling medium is assumed to be already set up. 
    Biomass formation is assumed to be already the objective. 
    
    model_id: name of the putput column (if None, 'output' will be used)
    biosynth: minimum fraction of the maximal biomass to guarantee.
    """    

    
    
    df = [] # list of dict to be converted in pnd dataframe
    if model_id == None: model_id = 'output'
    
    
    # get the max biomass
    res_before, obj_value_before, status_before = fba_no_warnings(model)
    
    
    for m in model.metabolites:
        if m.id.endswith('_c') == False:
            continue   # only interested in cytosolic metabolites

        with model:  # reversible changes 
            
            model.reactions.get_by_id(get_objectives(model)[0]).lower_bound = biosynth * obj_value_before
            
            binary, obj_value, status = can_synth(model, m.id)

            if status == 'optimal' and obj_value > 0.001:  
                able = 1 
            else:
                able = 0

            # save results in a future pnd DataFrame:
            df.append({'mid': m.id, model_id: able})
    
    
    df = pnd.DataFrame.from_records(df)
    df = df.set_index('mid', drop=True, verify_integrity=True)
    return df
    


def task_biosynth(accession, args):
    
    
    # retrive the arguments:
    outdir = args['outdir']
    skipgf = args['skipgf']
    biosynth = args['biosynth']
    
    
    # read json/sbml file:
    if not skipgf:
        ss_model = cobra.io.load_json_model(outdir + f'strain_models_gf/{accession}.json')
    else:  # user asked to skip the strain-specific gapfilling step
        ss_model = cobra.io.load_json_model(outdir + f'strain_models/{accession}.json')
    
    # perform the simulations: 
    df_results = biosynth_simulation(ss_model, model_id=accession, biosynth=biosynth)
    df_results = df_results.T.reset_index(drop=False).rename(columns={'index': 'accession'})
        
    # it has just 1 row:
    return [df_results.iloc[0].to_dict()]



def strain_biosynth_tests(logger, outdir, cores, panmodel, pam, skipgf, biosynth):
    
    
    # log some messages
    logger.info("Testing strain-specific biosynthesis of metabolites...")
    if biosynth != 0: 
        logger.info(f"Biomass will be fixed at {biosynth} fraction of its maximum.")
        
        
    # get the header for the rusults table:
    header = []
    for m in panmodel.metabolites: 
        if m.id.endswith('_c')==False:
            continue   # interested only in cytosolyc metabolites
        header.append(m.id)

    
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
            itertools.repeat(task_biosynth), 
            itertools.repeat({'outdir': outdir, 'skipgf': skipgf, 'biosynth': biosynth}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # save the table in the sae format of 'rpam':
    biosynth_pam = all_df_combined
    biosynth_pam = biosynth_pam.T
    biosynth_pam.to_csv(outdir + 'biosynth.csv')
    
    return 0