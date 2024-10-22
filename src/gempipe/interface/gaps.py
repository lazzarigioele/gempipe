import multiprocessing
import time
import os
from importlib import resources


from IPython.display import display
import pandas as pnd
import cobra
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.util.solver import linear_reaction_coefficients



__GAPSCACHE__ = None



def initialize(outdir):
    """Initialize the gempipe.curate API with the outputs coming from ``gempipe recon``.
    
    This function locates the draft pan-GSMM, PAM, and functional annotation table inside the `gempipe recon`` output folder (``-o``/``--outdir``).
    
    Args:
        outdir (str): path to the main output folder of ``gempipe recon`` (``-o``/``--outdir``). 
    
    Returns:
        cobra.Model: draft pan-GSMM to start the manual curation. 
    """
    
    
    # check the existance of the needed files: 
    if outdir.endswith('/')==False:
        outdir = outdir + '/'
    if not os.path.exists(outdir):
        print(f"ERROR: the specified path doesn't exists ({outdir}).")
        return
    panmodel_path = outdir + 'draft_panmodel.json'
    if not os.path.exists(panmodel_path):
        print(f"ERROR: cannot find a draft panmodel at the specified path ({panmodel_path}).")
        return
    pam_path = outdir + 'pam.csv'
    if not os.path.exists(pam_path):
        print(f"ERROR: cannot find a PAM at the specified path ({pam_path}).")
        return
    annot_path = outdir + 'annotation.csv'
    if not os.path.exists(annot_path):
        print(f"ERROR: cannot find a functional annotation table at the specified path ({annot_path}).")
        return
    report_path = outdir + 'report.csv'
    if not os.path.exists(report_path):
        print(f"ERROR: cannot find a report table at the specified path ({report_path}).")
        return
    
    
    # initilize or re-initialize the cache:
    global __GAPSCACHE__
    __GAPSCACHE__ = {}
    print(f"Loading PAM ({pam_path})...")
    __GAPSCACHE__['pam'] = pnd.read_csv(pam_path, index_col=0)
    print(f"Loading functional annotation table ({annot_path})...")
    __GAPSCACHE__['annot'] = pnd.read_csv(annot_path, index_col=0)
    print(f"Loading report table ({report_path})...")
    __GAPSCACHE__['report'] = pnd.read_csv(report_path, index_col=0)
    
    
    print(f"Loading draft pan-GSMM ({panmodel_path})...")
    return cobra.io.load_json_model(panmodel_path)



def get_objectives(model):
    """Get the IDs of the current objective reactions. 
    
    Args:
        model (cobra.Model): target model.
        
    Returns:
        list: IDs of the reactions set as objective.

    """
    
    objs = list(linear_reaction_coefficients(model).keys())
    obj_ids = [obj.id for obj in objs]
    return obj_ids

        
        
def get_solver(model):
    """Get the ID of the solver associated to the model.
    
    Args:
        model (cobra.Model): target model.
        
    Returns:
        str: ID of the solver (for example: ``glpk_exact``).

    """
    
    solver = str(type(model.solver))
    solver = solver.replace('<class ', '')
    solver = solver.replace("'optlang.", '')
    solver = solver.replace("_interface.Model'>", '')
    
    return solver



def remove_rids(model, rids, inverse=False):
    """Remove reactions from the model given a list of reaction IDs.
    
    Args:
        model (cobra.Model): target model.
        rids (list): reaction IDs.
        inverse (bool): if ``True``, reactions IDs contained in `rids` will be the ones to keep and not to remove.

    """
    
    to_delete = []
    for r in model.reactions: 
        if not inverse:
            if r.id in rids:
                to_delete.append(r)
        else:
            if r.id not in rids:
                to_delete.append(r)
    model.remove_reactions(to_delete)


        
def perform_gapfilling(model, universe, mid=None, slim=None, minflux=1.0, exr=False, nsol=3, penalties=None, verbose=True, timeout=None, logger=None, boost=False): 
    """Propose gap-filling solutions for the specified objective. 
    
    It's possible to gap-fill also for the biosynthesis of a specific metabolite.
    
    Args:
        model (cobra.Model): target model to gap-fill.
        universe (cobra.Model): model from which to take new reactions.
        mid (str): gap-fill for the biosynthesis of a specific metabolite having ID `mid`. Will be ignored if ``None``.
        slim (str): try to reduce the complexity of the universe, considering only its reactions carrying non-0 flux. Can be ``FBA`` or ``FVA``. Will be ignored if ``None``.
        minflux (float): minimal flux to grant through the objective reaction.
        nsol (int): number of alternative solutions. 
        exr (bool): whether to allow the opening of new EX_change reactions.
        penalties (dict): dictionary keyed by reaction ID, containing reaction-specific penalties to apply during gap-filling.
        verbose (bool): if False, just return the lisr of reaction IDs without printing any further information.
        timeout (int): max seconds to wait for the gapfilling step. If ``None``, gap-filling won't be temporized. 
        logger (logging.Logger): write exception on a logger instad of using print().
        boost (bool): if ``True``, consider the current nutritive sources as unlimited (for example, if LB of ``EX_fe3_e`` is -0.0075, it will be raised to -1000).
        
    Returns:
        list: IDs of reactions proposed during the 1st solution.
    """
    
    
    def temporized_task(model, universe, minflux, exr, nsol, penalties):

        # launch cobrapy gapfilling: 
        gapfiller = GapFiller(model, universe, 
            lower_bound = minflux,
            demand_reactions = False, 
            exchange_reactions = exr, 
            integer_threshold = 1e-10,  # default: 1e-6.
            penalties = penalties,
        )
        try: solutions = gapfiller.fill(iterations=nsol)
        except Exception as e: 
            return e
        
        return solutions
            
            
    # temporary changes (objective, solver, demands edits are temporary)
    with model, universe: 

        """
        # set new solver if needed (cannot be 'glpk_exact').
        if get_solver(model) != solver: model.solver = solver
        if get_solver(universe) != solver: universe.solver = solver
        """

        # if focusing on a particular biosynthesis
        if mid != None:
            model.objective = add_demand(model, mid)
            universe.objective = add_demand(universe, mid)
            
            
        # if requested, boost the actual sources:
        for r in model.reactions: 
            if len(r.metabolites)==1 and list(r.metabolites)[0].id.endswith('_e'): # if exchange reaction
                if r.lower_bound < 0: 
                    r.lower_bound = -1000


        # if requested, try to reduce the universe complexity:
        if slim == 'FBA':
            fluxes = universe.optimize().fluxes
            rids_to_keep = fluxes[fluxes != 0].index
            remove_rids(universe, rids_to_keep, inverse=True)
        elif slim == 'FVA':
            fluxes = cobra.flux_analysis.flux_variability_analysis(universe, fraction_of_optimum=0.01, loopless=False)
            rids_to_keep = fluxes[(fluxes['minimum']!=0) | (fluxes['maximum']!=0)].index
            remove_rids(universe, rids_to_keep, inverse=True)
            
            
        # "AssertionError: daemonic processes are not allowed to have children"
        # could be raise as child processes created by multiprocessing.Pool are demoniac by deafault,
        # and demoniac processes cannot create their own child processes like the temporized gap-filler. 
        if timeout != None: 
            
            
            # run the gapfilling algorithm inside a timer: 
            # first define the worker for the separate process
            def worker(results_channel):
                result = temporized_task(model, universe, minflux, exr, nsol, penalties)
                results_channel.put(result)


            # then start the separate process
            results_channel = multiprocessing.Queue()
            process = multiprocessing.Process(target=worker, args=(results_channel,))
            process.start()  # start the process
            process.join(timeout) # wait for the process to finish or timeout
        
        
            # if still running after the timeout:
            if process.is_alive():  
                process.terminate()
                e = "Gap-filling timed out, consider to increase the 'timeout' attribute."
                if verbose:  # avoid the error stack trace
                    print("ERROR: cobrapy:", e) if logger==None else logger.debug(f'cobrapy: {e}')  
                    #raise TimeoutError(result)  # show the error trace
                return None
        
        
            else: # the worker finished: retrieve the result or exception from the queue:
                result = results_channel.get()
                if isinstance(result, Exception):  
                    if verbose: 
                        print("ERROR: cobrapy:", result) if logger==None else logger.debug(f'cobrapy: {result}')  
                        #raise Exception(result)  # show the error trace
                    return None  # exit the function
                else:  # a true result
                    solutions = result
                    
                    
        else:  # do not temporize the gap-filling:
            result = temporized_task(model, universe, minflux, exr, nsol, penalties)
            if isinstance(result, Exception):  
                if verbose: 
                    print("ERROR: cobrapy:", result) if logger==None else logger.debug(f'cobrapy: {result}')
                return None
            else:  # a true result
                solutions = result
                

        # iterate the solutions:
        first_sol_rids = []  # rids proposed during the 1st solution
        for i, solution in enumerate(solutions):
            if verbose and logger==None: print(f'Solution {i+1}. Reactions to add: {len(solution)}.')


            # iterate the reactions: 
            counter = 0
            for r in solution: 
                counter += 1
                # Note: this 'r' is not linked to any model. 
                # Indeed, calling r.model, None will be returned. 
                if verbose and logger==None: print(f'{counter} {r.id} {r.name}')

                # populate results with IDs from first solution:
                if i == 0: first_sol_rids.append(r.id)

            # separate solutions with a new line:
            if i+1 != len(solutions): 
                if verbose and logger==None: print()


        return first_sol_rids
        
    

def get_universe(staining='neg'):
    """Return a CarveMe universe. 
    
    Args:
        staining (str): 'pos' or 'neg'.
        
    Returns: 
        cobra.Model: the selected universe.
    """
    
    # basically it's a wrapper of the recon function
    from gempipe.recon.networkrec import get_universe_template
    universe = get_universe_template(logger=None, staining=staining)
    
    return universe



def get_biolog_mappings():
    """Return the Biolog mappings internally used by gempipe.
    
    Plate information is taken from DuctApe (https://doi.org/10.1016/j.ygeno.2013.11.005).
        
    Returns: 
        pandas.DataFrame: the Biolog mappings.
    """
    
    with resources.path('gempipe.assets', 'biolog_mappings.csv' ) as asset_path: 
        biolog_mappings = pnd.read_csv(asset_path, index_col=0)
    
    return biolog_mappings



def add_demand(model, mid):
    """Create a demand reaction, useful for debugging models.
    
    Args:
        model (cobra.Model): target model.
        mid (str): metabolite ID (compartment included) for which to create the demand.
        
    Returns:
        str: demand reaction ID.
    """
    
    rid = f"demand_{mid}"
    newr = cobra.Reaction(rid)
    model.add_reactions([newr])
    model.reactions.get_by_id(rid).reaction = f"{mid} -->"
    model.reactions.get_by_id(rid).bounds = (0, 1000)
    
    return rid



def can_synth(model, mid):
    """Check if the model can synthesize a given metabolite.
    
    Args:
        model (cobra.Model): target model.
        mid (str): metabolite ID (compartment included) for which to check the synthesis.
    
    Returns:
        (bool, float, str):
        
            `[0]` ``True`` if `mid` can be synthesized (``optimal`` status and positive flux).
        
            `[1]` maximal theoretical synthesis flux.
            
            `[2]` status of the optimizer.
    """
    
    # changes are temporary: demand is not added, objective is not changed.
    with model: 

        rid = add_demand(model, mid)

        # set the objective to this demand reaction:
        model.objective = rid

        # perform FBA: 
        res = model.optimize()
        value = res.objective_value
        status = res.status
        response = True if (value > 0 and status == 'optimal') else False
        
        return response, round(value, 2), status
    
    
    
def check_reactants(model, rid, verbose=True):
    """Check which reactant of a given reaction cannot be synthesized.
    
    Args:
        model (cobra.Model): target model.
        rid (str): reaction ID for which to check the synthesis of the reactants.
        verbose (bool): if False, just return the list of metabolite IDs without printing any further information.

    
    Returns:
        list: IDs of blocked reactants.
    """
    
    # changes are temporary
    with model: 
        counter = 0

        
        # get reactants and products
        reacs = [m for m in model.reactions.get_by_id(rid).reactants]
        prods = [m for m in model.reactions.get_by_id(rid).products]

        
        # iterate through the reactants: 
        mid_blocked = []
        for i, m in enumerate(reacs):
            
            # check if it can be synthesized:
            response, flux, status = can_synth(model, mid=m.id)
            
            if response==False: 
                counter += 1
                if verbose: print(f'{counter} : {flux} : {status} : {m.id} : {m.name}')
                mid_blocked.append(m.id)

                
        return mid_blocked
    
    
    
def sensitivity_analysis(model, scaled=False, top=3, mid=None):
    """Perform a sensitivity analysis (or reduced costs analysis) focused on the EX_change reaction.
    
    It is based on the current model's objective. The returned dictionary is sorted from most negative to most positive values.
    
    Args:
        model (cobra.Model): target model.
        scaled (bool): whether to scale to the current objective value.
        top (int): get just the first and last `top` EX_change reactions. If ``None``, all EX_change reactions will be returned.
        mid (str): instead of optimizing for the current objective reaction, do the analysis on the biosynthesis of a specific metabolite having ID `mid`. If `None` it will be ignored.
    
    Returns:
        dict: reduced costs keyd by EX_change reaction ID. 
    """
    
    
    # temporary chenges:
    with model:
        
      
        # focus on a specific metbaolite
        if mid != None: 
            model.objective = add_demand(model, mid)


        res = model.optimize()
        obj = res.objective_value
        flx = res.fluxes.to_dict()
        rcs = res.reduced_costs.to_dict()


        # manage 0 flux exception:
        if obj == 0 and scaled == True:
            raise Exception("Cannot scale reduced costs id the objective value is 0")


        # get the reduced costs of the EXR only:
        rcsex = {} 
        for key in rcs:
            if key.startswith("EX_"):
                if not scaled : rcsex[key] = rcs[key]
                else: rcsex[key] = rcs[key] * flx[key] / obj


        # get the most impactful (lowest and highest)
        rcsex = sorted(rcsex.items(), key=lambda x: x[1], reverse=True)
        rcsex = {i[0]: i[1] for i in rcsex}  # convert list to dictionary
    
        
        # get only the top N and bottom N exchanges
        if top != None: 
            rcsex_filt = {}
            for i in range(top):
                rcsex_filt[ list(rcsex.keys())[i] ] = list(rcsex.values())[i]
            for i in range(top):
                rcsex_filt[ list(rcsex.keys())[-top +i] ] = list(rcsex.values())[-top +i]
            rcsex = rcsex_filt
        
        return rcsex
    
    
    
def query_pam(name=[], ko=[], kr=[], km=[], kt=[], ec=[], des=[], pfam=[], annot=False, model=None, rid=None):
    """Show clusters in the context of a PAM.
    
    Clusters can be selected based on thier functional annotation.
    
    Args: 
        name (list): preferred names to search for, eg ['fabB', 'fabG'].
        ko (list): KEGG Orthologs (KOs) to search for, eg ['K00647', 'K00059'].
        kr (list): KEGG Reaction IDs to search for, eg ['R03460'].
        km (list): KEGG Module IDs to search for, eg ['M00846'].
        kt (list): KEGG Transport IDs to search for, eg ['3.A.3.2'].
        ec (list): EC codes to search for, eg ['2.5.1.19'].
        des (list): gene function descriptions to search for.
        pfam (list): PFAM domains to search for, eg ['Pyridox_ox_2'].
        annot (bool): if ``True``, return the functional annotation table instead of the PAM. 
        model (cobra.Model): if not ``None``, create a new column `modeled` as first, revealing the presence of each cluster in the model.
        rid (list): show clusters involved in the specified reactions (provided their IDs). Requires ``model``.
    
    Returns:
        pandas.DataFrame: filtered PAM or annotation table. 
    """
    
    
    # get the needed file from the cache:
    global __GAPSCACHE__
    if __GAPSCACHE__ == None:
        print("ERROR: you first need to execute gempipe.initialize().")
        return
    annotation = __GAPSCACHE__['annot']
    pam = __GAPSCACHE__['pam']
    report = __GAPSCACHE__['report']

    
    # create a copy to filter: 
    annotation_filter = annotation.copy()
    cluster_set = set(annotation_filter.index.to_list())
    
    
    # filter based on involved clusters: 
    if model != None: 
        if rid != None: 
            if type(rid)==str: rid = [rid]
            good_clusters = []
            for i in rid: 
                r = model.reactions.get_by_id(i)
                good_clusters = good_clusters + [g.id for g in r.genes if g.id in annotation.index]
            cluster_set = cluster_set.intersection(set(good_clusters))
                
    
    # filter for kegg orthologs
    if ko != []:
        if type(ko)==str: ko = [ko]
        good_clusters = []
        for i in ko:
            good_clusters = good_clusters + list(annotation[annotation['KEGG_ko'].str.contains(f'ko:{i}')].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
       
    # filter for kegg reactions
    if kr != []:
        if type(kr)==str: kr = [kr]
        good_clusters = []
        for i in kr:
            good_clusters = good_clusters + list(annotation[annotation['KEGG_Reaction'].str.contains(f'{i}')].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
        
    # filter for KEGG Modules
    if km != []:
        if type(km)==str: km = [km]
        good_clusters = []
        for i in km:
            good_clusters = good_clusters + list(annotation[annotation['KEGG_Module'].str.contains(f'{i}')].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
        
    # filter for KEGG Transport
    if kt != []:
        if type(kt)==str: kt = [kt]
        good_clusters = []
        for i in kt:
            good_clusters = good_clusters + list(annotation[annotation['KEGG_TC'].str.contains(f'{i}', case=False)].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
    
    
    # filter for gene symbols
    if name != []:
        if type(name)==str: name = [name]
        good_clusters = []
        for i in name:
            good_clusters = good_clusters + list(annotation[annotation['Preferred_name'].str.contains(f'{i.lower()}', case=False)].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
        
    # filter for EC:
    if ec != []:
        if type(ec)==str: ec = [ec]
        good_clusters = []
        for i in ec:
            good_clusters = good_clusters + list(annotation[annotation['EC'].str.contains(f'{i}', case=False)].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
        
    # filter for pfam domains
    if pfam != []:
        if type(pfam)==str: pfam = [pfam]
        good_clusters = []
        for i in pfam:
            good_clusters = good_clusters + list(annotation[annotation['PFAMs'].str.contains(f'{i}', case=False)].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
        
    # filter for function description
    if des != []:
        if type(des)==str: des = [des]
        good_clusters = []
        for i in des:
            good_clusters = good_clusters + list(annotation[annotation['Description'].str.contains(f'{i}', case=False)].index)
        cluster_set = cluster_set.intersection(set(good_clusters))
        
    
    # get tabular results (PAM or annotation)
    annotation_filter = annotation_filter.loc[list(cluster_set), ]
    if annot: results = annotation_filter
    else: results = pam.loc[[i for i in annotation_filter.index if i in pam.index], ]
    
    
    # rename columns using the species information:
    if not annot:
        old_to_new = {}
        for column in results.columns:
            match = report[report['accession']==column].iloc[0]
            old_to_new[column] = f"{match['species']} {match['accession']}"
        results = results.rename(columns=old_to_new)  # rename columns.
        results = results.reindex(sorted(results.columns), axis=1)  # alphabetical order.
    
    
    # mark clusters that are already modeled:
    if model != None: 
        results_columns = list(results.columns)
        results['modeled'] = 'False'
        for cluster in results.index:
            if cluster in [g.id for g in model.genes]:
                results.loc[cluster, 'modeled'] = ', '.join([r.id for r in model.genes.get_by_id(cluster).reactions])
        results = results[['modeled'] + results_columns]  # reorder columns
    
    
    return results



def import_from_universe(model, universe, rid, bounds=None, gpr=None):
    """Insert a new reaction taken from a universe model.

    Args:
        model (cobra.Model): target model to expand with new reactions.
        universe (cobra.Model): universe model, source of new reactions.
        rid (str): Id of the reaction to transfer. 
        bounds (tuple): bounds to apply to the inserted reaction, eg (0, 1000). If ``None``, bounds from universe will be retained.
        gpr (str): GPR to associate to the inserted reaction. If ``None``, no GPR will be associated.
    
    """
    
    r = universe.reactions.get_by_id(rid)
    model.add_reactions([r])
    
    if bounds != None:
        model.reactions.get_by_id(rid).bounds = bounds
        
    if gpr != None:
        model.reactions.get_by_id(rid).gene_reaction_rule = gpr
        model.reactions.get_by_id(rid).update_genes_from_gpr()
        
        if gpr == '':
            # be shure we are not importing panmodel genes:
            to_remove = []
            involved_gids = [g.id for g in model.genes]
            for g in model.genes: 
                if g.id in involved_gids and len(g.reactions)==0:
                    to_remove.append(g)
            cobra.manipulation.delete.remove_genes(model, to_remove, remove_reactions=True)



def ss_preview(panmodel, accession):
    """Create a preview of a strain-specific model given its accession and the starting draft pan-GSMM.
    
    Args:
        panmodel (cobra.Model): draft pan-GSMM (ideally undergoing manual curation) from which to derive the strain-specific model.
        accession (str): accession of the strain for which to derive the strain-specific model.
        
    Returns:
        cobra.Model: strain-specific metabolic model.
        
    """
    
    # import here to avoid the "circular import" error:
    from gempipe.derive.strain import subtract_clusters
    
    
    # get the needed file from the cache:
    global __GAPSCACHE__
    if __GAPSCACHE__ == None:
        print("ERROR: you first need to execute gempipe.initialize().")
        return
    annotation = __GAPSCACHE__['annot']
    pam = __GAPSCACHE__['pam']
    report = __GAPSCACHE__['report']
    
    
    # create strain specific model
    ss_model = panmodel.copy()  
    response = subtract_clusters(ss_model, panmodel, pam, accession)
    
    
    return ss_model



def biolog_preview(model, starting_C=None, starting_N=None, starting_P=None, starting_S=None, seed=False):
    """Get a preview of of the Biolog simulations for the provided model.
    
    Args:
        model (cobra.Model): target model.
        starting_C (str): starting C source; if `None`, will be assigned to `EX_glc__D_e` or `EX_cpd00027_e0` depending on the selected database.
        starting_N (str): starting N source; if `None`, will be assigned to `EX_nh4_e` or `EX_cpd00013_e0` depending on the selected database.
        starting_P (str): starting P source; if `None`, will be assigned to `EX_pi_e` or `EX_cpd00009_e0` depending on the selected database.
        starting_S (str): starting S source; if `None`, will be assigned to `EX_so4_e` or `EX_cpd00048_e0` depending on the selected database.
        seed (bool): if `True`, use the SEED notation instead of the BiGG notation.
        
    Returns:
        pandas.DataFrame: Biolog simulations.
        
    """
    
    # import here to avoid the "circular import" error:
    from gempipe.derive.biolog import biolog_simulation
    
    
    if starting_C==None: 
        if seed==True: starting_C = 'EX_cpd00027_e0'
        else: starting_C = 'EX_glc__D_e'
    if starting_N==None: 
        if seed==True: starting_N = 'EX_cpd00013_e0'
        else: starting_N = 'EX_nh4_e'
    if starting_P==None: 
        if seed==True: starting_P = 'EX_cpd00009_e0'
        else: starting_P = 'EX_pi_e'
    if starting_S==None: 
        if seed==True: starting_S = 'EX_cpd00048_e0'
        else: starting_S = 'EX_so4_e'
    
    # do the simulations
    biolog_mappings = get_biolog_mappings()
    results_df = biolog_simulation(model, biolog_mappings, seed, starting_C, starting_N, starting_P, starting_S)
    
    
    return results_df



def add_reaction(model, rid, rstring, bounds=None, name='', gpr=None, ):
    """Add a new reaction to the model.
    
    Useful for example when the universe needs to be expanded.
    
    Args:
        model (cobra.Model): target model to expand with the new reaction.
        rid (str): ID for the new reaction.
        rstring (str): string describing the reaction. 
        bounds (set): lower and upper bound for the new reaction. 
        name (str): a name for the new reaction. 
        gpr (str): string describing the GPR.
        
    """
    
    r = cobra.Reaction(rid)
    r.name = name
    model.add_reactions([r])
    r = model.reactions.get_by_id(rid)
    r.build_reaction_from_string(rstring)
    if bounds != None: 
        r.bounds = bounds
    if gpr != None: 
        r.gene_reaction_rule = gpr
        r.update_genes_from_gpr()
        
        

def add_metabolite(model, mid, formula=None, charge=None, name=''):
    """Add a new metabolite to the model.
    
    Useful for example when the universe needs to be expanded.
    
    Args:
        model (cobra.Model): target model to expand with the new metabolite.
        mid (str): ID for the new metabolite.
        formula (str): string the chemical formula.
        charge (int): charge of the metabolite. 
        name (str): a name for the new metabolite. 
        
    """
    
    m = cobra.Metabolite(mid)
    m.name = name
    model.add_metabolites([m])
    m = model.metabolites.get_by_id(mid)
    m.formula = formula
    m.charge = charge
    m.compartment = mid.rsplit('_')[-1]
        


def search_similar(panmodel, rid, field='ko', unmod=False, species=None, showgpr=False, forceshow=False):
    """Search the PAM for gene clusters having functional annotation similar to that of gene clusters involved in the specified reaction of the pan-GSMM.
    
    First, the reaction `rid` is extracted from the `panmodel`. 
    Then, gene clusters are extracted from its GPR.
    The PAM is finally searched for all the gene clusters having the same functional annotation in one of the following fields: 'ko', 'pfam', 'kt', and 'ec'.
    
    Args:
        panmodel (cobra.Model): a pan-GSMM having the reaction `rid`.
        rid (str): ID of the reaction using the gene clusters on which to focus the search.
        field (str): functional annotation field to be used to extract similar gene clusters.
        unmod (bool): if ``True``, show only gene clusters that are not yet part of the `panmodel` (unmodeled).
        species (str): show only columns (genomes) having this substring in thier name (useful to filter for particular species or strains).
        showgpr (bool): if ``True``, show the pan-GSMM GPR of the selected reaction `rid`.
        forceshow (bool): if ``True``, display the results table with no limits on the number of displayable rows and columns. 
    
    Returns:
        pandas.DataFrame: PAM-derived results table, mixing presence/absence of gene clusters together with their functional annotation.

    """
    
    # check the selected feat: 
    available_fields = ['ko', 'pfam', 'kt', 'ec']
    if field not in available_fields: 
        print(f"field not recognized. Please use one of {available_fields}.")
        return
    results = None
    
    
    # show GPR in panmodel
    if showgpr: 
        print(f'GPR for {rid}:', panmodel.reactions.get_by_id(rid).gene_reaction_rule)
    involved_genes = [g.id for g in panmodel.reactions.get_by_id(rid).genes]
    
    
    def expand_terms(terms):
        # expand terms as they could be multiple (eg: ['K01620', 'K01620,K20801'])
        terms_exp = [] 
        for i in terms: 
            terms_exp = terms_exp + i.split(',')
        terms_exp = list(set(terms_exp))
        return terms_exp
    
    
    if field == 'pfam': 
        terms = list(set(query_pam(annot=True).loc[involved_genes, ]['PFAMs'].to_list()))
        if terms == ['-']: return None
        terms = expand_terms(terms)
        presence_absence = query_pam(pfam=terms, model=panmodel)
        func_annot = query_pam(pfam=terms, annot=True)[['Description', 'KEGG_ko', 'PFAMs']]
        
        
    if field == 'ko': 
        terms = list(set(query_pam(annot=True).loc[involved_genes, ]['KEGG_ko'].to_list()))
        if terms == ['-']: return None
        terms = [i.replace('ko:', '') for i in terms]
        terms = expand_terms(terms)
        presence_absence = query_pam(ko=terms, model=panmodel)
        func_annot = query_pam(ko=terms, annot=True)[['Description', 'KEGG_ko', 'PFAMs']]
        
        
    if field == 'kt': 
        terms = list(set(query_pam(annot=True).loc[involved_genes, ]['KEGG_TC'].to_list()))
        if terms == ['-']: return None
        terms = expand_terms(terms)
        presence_absence = query_pam(kt=terms, model=panmodel)
        func_annot = query_pam(kt=terms, annot=True)[['Description', 'KEGG_TC', 'PFAMs']]
        
        
    if field == 'ec': 
        terms = list(set(query_pam(annot=True).loc[involved_genes, ]['EC'].to_list()))
        if terms == ['-']: return None
        terms = expand_terms(terms)
        presence_absence = query_pam(ec=terms, model=panmodel)
        func_annot = query_pam(ec=terms, annot=True)[['Description', 'EC', 'KEGG_ko', 'PFAMs']]
            
    
    # filter for desired columns: 
    if species != None: 
        presence_absence = presence_absence[['modeled']+[i for i in presence_absence.columns if species in i]]

    
    # merge presence/absence table and functional annootation table:
    results = pnd.concat([presence_absence, func_annot], axis=1)
    
    
    # exclude gene clusters already modeled: 
    if unmod: 
        results = results[results['modeled']=='False']
    
    
    # compute percentage of genomes with no representative (among this selection of gene clusters)
    empty_columns = results.columns[results.isnull().all()]
    structure_cols = ['modeled', 'Description', 'EC', 'KEGG_ko', 'PFAMs', 'KEGG_TC']  # cols not related to input genomes
    empty_columns = set(empty_columns) - set(structure_cols)
    all_columns = set(list(results.columns)) - set(structure_cols)
    print(f'Non-empty columns: {len(all_columns)-len(empty_columns)}/{len(all_columns)} ({round((len(all_columns)-len(empty_columns))/len(all_columns)*100,1)}%); {field} terms considered: {str(terms)}.')

    
    if forceshow:   # temporary allow pandas/jupyter to display with no row/col limits: 
        with pnd.option_context('display.max_rows', None, 'display.max_columns', None):
            display(results)
    
    return results