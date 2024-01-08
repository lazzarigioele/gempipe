import multiprocessing
import time


import cobra
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.util.solver import linear_reaction_coefficients



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


        
def perform_gapfilling(model, universe, minflux=1.0, nsol=3, add=False, solver='glpk', exr=False, timeout=120): 
    """Propose gap-filling solutions for the specified objective. 
    
    Args:
        model (cobra.Model): target model to gap-fill.
        universe (cobra.Model): model from which to take new reactions.
        minflux (float): minimal flux to grant through the objective reaction.
        nsol (int): number of alternative solutions. 
        add (bool): if ``True``, suggested reactions are added into `model`, otherwise they are simply listed.
        solver (str): solver to use (usually ``glpk`` or ``cplex``).
        exr (bool): whether to allow the opening of new EX_change reactions.
        timeout (int): maximum amount of seconds to wait for this function to finish.
        
    """
    
    
    def inner_gafiller(model, universe, minflux, nsol, add, solver, exr):
        
        
        # get the original solvers
        ori_solver_model = get_solver(model)
        ori_solver_universe = get_solver(model)


        # set new solver if needed:
        # solver cannot be 'glpk_exact'.
        if ori_solver_model != solver: model.solver = solver
        if ori_solver_universe != solver: universe.solver = solver


        # compute reactions to add.
        gapfiller = GapFiller(model, universe, 
            lower_bound = minflux,
            demand_reactions = False,
            exchange_reactions = exr,
            integer_threshold = 1e-9,  
        )
        
        solutions = gapfiller.fill(iterations=nsol)
        
        
        # iterate the solutions:
        for i, solution in enumerate(solutions):
            print(f'Solution {i+1}. Reactions to add: {len(solution)}.')


            # iterate the reactions: 
            counter = 0
            for reaction in solution: 
                counter += 1
                print(counter, reaction.id, reaction.name)

                # Add the reaction if requested.
                # Note: this reaction is not linked to any model. In fact, (1) if you call
                # reaction.model, None will be returned. Instead, if we do 
                # universe.reactions.get_by_id(reaction.id).model, we obtain 'pan_lpl'
                if add: model.add_reactions([reaction])

            
        # restore the original solver if needed:
        if ori_solver_model != solver: model.solver = ori_solver_model
        if ori_solver_universe != solver: universe.solver = ori_solver_universe
        
        
        return None
        
        
    
    # run the gapfilling algorithm inside a timer: 
    def worker(results_channel):
        try:
            result = inner_gafiller(model, universe, minflux, nsol, add, solver, exr)
            results_channel.put(result)
        except Exception as exception:  # eventual exception raised by the gapfilling algorithm itself
            results_channel.put(exception)  # put exceptions inside the same channel of results.
    results_channel = multiprocessing.Queue()
    process = multiprocessing.Process(target=worker, args=(results_channel,))
    process.start()  # start the process
    process.join(timeout) # wait for the process to finish or timeout
    if process.is_alive():  # if still running after the timeout
        process.terminate()
        process.join()
        raise TimeoutError("perform_gapfilling() timed out: please consider the perform_gapfilling_slim() alternative.")
    else: # retrieve the result or exception from the queue
        result = results_channel.get()
        if isinstance(result, Exception):
            print("ERROR:", result)  # print to avoid the annoing error trace
            #raise Exception(result)  # show the error trace
            return  # exit the function
        else:  # a true result
            return result
        
    
    

    

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
        verbose (bool): whether to print blocked reactants.
    
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
        verbose_string = ''
        for m in reacs:
            
            # check if it can be synthesized:
            response, flux, status = can_synth(model, mid=m.id)
            
            if response==False: 
                counter += 1
                verbose_string = verbose_string + f'{counter} : {flux} : {status} : {m.id} : {m.name}\n'
                mid_blocked.append(m.id)

                
        verbose_string = verbose_string.rstrip()  # remove trailing endline
        if verbose: print(verbose_string)
        return mid_blocked