import re
import importlib.metadata


import cobra


from gempipe.interface.gaps import get_solver
from gempipe.interface.gaps import get_objectives


def close_boundaries(model):
    """Set all the EX_change reactions to (0, 0).
    
    Args:
        model (cobra.Model): target model.

    """
    
    for r in model.reactions:
        if len(r.metabolites)==1 and r.id.startswith('EX_'):
            r.bounds = (0, 0)



def verify_egc(model, mid, escher=False, threshold=1e-5, verbose=True): 
    """Test the presence of energy-generating cycles (EGCs). 
    
    Can also output a model for Escher, with just the reactions composing the cycle. 
    
    Args:
        model (cobra.Model): target model. Must be encoded with the BiGG notation.
        mid (str): metabolite ID for which the EGC must be checked. Warning: must be without compartment, so for example ``atp`` instead of ``atp_c``. 
        escher (bool): save a reduced ``cobra.Model`` in the current directory. To be loaded in Escher. 
        threshold (float): values below this treshold are considered as 0.
        verbose (bool): if True, show results of the test even if no EGC was detected.
        
    Returns:
        bool: ``True`` if an EGC is detected.
    """
    
    # changes as not permament: 
    with model: 
        
        # close all the exchange reactions: 
        close_boundaries(model)

                
        # create a dissipation reaction: 
        dissip = cobra.Reaction(f'__dissip__{mid}')
        model.add_reactions([dissip])
        dissip = model.reactions.get_by_id(f'__dissip__{mid}')
        
        
        # define the dissipation reaction:
        modeled_mids = [m.id for m in model.metabolites]
        if mid == 'atp':
            dissip_string = 'atp_c + h2o_c --> adp_c + pi_c + h_c'
        elif mid == 'ctp':
            dissip_string = 'ctp_c + h2o_c --> cdp_c + pi_c + h_c'
        elif mid == 'gtp':
            dissip_string = 'gtp_c + h2o_c --> gdp_c + pi_c + h_c'
        elif mid == 'utp':
            dissip_string = 'utp_c + h2o_c --> udp_c + pi_c + h_c'
        elif mid == 'itp':
            dissip_string = 'itp_c + h2o_c --> idp_c + pi_c + h_c'
        elif mid == 'nadh':
            dissip_string = 'nadh_c --> nad_c + h_c'
        elif mid == 'nadph':
            dissip_string = 'nadph_c --> nadp_c + h_c'
        elif mid == 'fadh2':
            dissip_string = 'fadh2_c --> fad_c + 2.0 h_c'
        elif mid == 'accoa':
            dissip_string = 'accoa_c + h2o_c --> ac_c + coa_c + h_c'
        elif mid == 'glu__L':
            if 'nh4_c' in modeled_mids :
                dissip_string = 'glu__L_c + h2o_c --> akg_c + nh4_c + 2.0 h_c'
            elif 'nh3_c' in modeled_mids :
                dissip_string = 'glu__L_c + h2o_c --> akg_c + nh3_c + 3.0 h_c'
            else:
                raise Exception("'nh4_c' or 'nh3_c' must be present in the model.")
        elif mid == 'q8h2':
            dissip_string = 'q8h2_c --> q8_c + 2.0 h_c'

        else: 
            raise Exception("Metabolite ID (mid) not recognized.") 
        dissip.build_reaction_from_string(dissip_string)
        
        
        # set the objective and optimize: 
        model.objective = f'__dissip__{mid}'
        res = model.optimize()
        
        
        # apply the threshold:
        obj_value = res.objective_value
        if abs(obj_value) < threshold:
            obj_value = 0
        
        
        # log some messages
        if obj_value != 0 or verbose:
            print(dissip.reaction)
            print(obj_value , ':', res.status )


            # get suspect !=0 fluxes (if any)
            if obj_value != 0: 
                fluxes = res.fluxes
                print()  # skip a line befor printing the EGC members

                # get interesting fluxes (0.001 tries to take into account the approximation in glpk and cplex solvers)
                fluxes_interesting = fluxes[((fluxes > 0.001) | (fluxes < -0.001)) & (fluxes.index != f'__dissip__{mid}')]
                print(fluxes_interesting.to_string())


                # create a model for escher
                if escher:  
                    model_copy = model.copy()
                    all_rids = [r.id for r in model_copy.reactions]
                    to_delete = set(all_rids) - set(fluxes_interesting.index)
                    model_copy.remove_reactions(to_delete)
                    cobra.io.save_json_model(model_copy, f'__dissip__{mid}' + '.json')
                    print(f'__dissip__{mid}' + '.json', "saved in current directory.")
                    
        if obj_value != 0: 
            return True
        else: return False



def verify_egc_all(model):
    """Quickly check the presence of EGCs for the main metabolites.
    
    Internally calls `verify_egc()` over a list of metabolite IDs.
    
    Args:
        model (cobra.Model): target model.
    """
    
    mids_to_check = ['atp','ctp','gtp','utp','itp','nadh','nadph','fadh2','accoa','glu__L','q8h2']
    all_results = []
    for mid in mids_to_check:
        result = verify_egc(model, mid, verbose=False)
        all_results.append(~result)
    if all(all_results):
        print("No energy-generating cycles (EGCs) found.")
        
        
                

    
def check_sink_demand(model, verbose=True):
    """Check presence of sink and demand reactions.
    
    Here they are simply defined as reactions involving just 1 metabolite, with an ID not starting with ``EX_``.
    
    Args:
        model (cobra.Model): target model.
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of sink/demand reactions found.
    """
    
    found_rids = []
    cnt = 0
    for r in model.reactions:
        if len(r.metabolites) == 1: 
            if r.id.startswith("EX_") == False: 
                cnt += 1
                if verbose: print(cnt, ':', r.id, ':', r.reaction, ':', r.bounds)
                found_rids.append(r.id)
    
    
    if found_rids == []:
        if verbose: print("No sink/demand reactions found.")
    return found_rids
   

                
def check_exr_notation(model, verbose=True): 
    """Check that every EX_change reaction ID begins with ``EX_``.
    
    Here EX_change reactions are defined as those reactions having just 1 metabolite involved, included in the extracellular compartment.
    
    Args:
        model (cobra.Model): target model. 
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of EX_change reactions with bad ID.
    """

    found_rids = []
    cnt = 0
    for r in model.reactions:
        if len(r.metabolites) == 1: 
            if list(r.metabolites)[0].id.rsplit('_', 1)[1] == 'e':  # extracellular compartment
                if r.id.startswith("EX_") == False: 
                    cnt += 1
                    if verbose: print(cnt, ':', r.id, ':', r.reaction)
                    found_rids.append(r.id)
    
    
    if found_rids == []:
        if verbose: print("No EX_change reaction with bad ID found.")
    return found_rids

 
    
def remove_EX_annots(model): 
    """Remove all annotations from EX_change reactions.
    
    Args:
        model (cobra.Model): target model. 
    """
    
    cnt = 0
    for r in model.reactions: 
        if len(r.metabolites) == 1 and r.id.startswith("EX_"): 
            if r.annotation != {}: 
                r.annotation = {}

    
    
def check_missing_charges(model, verbose=True):
    """Check if all metabolites have a charge attribute.
    
    Args:
        model (cobra.Model): target model. 
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of matabolites missing the charge attribute.
    
    """ 
    
    found_mids = []
    cnt = 0
    for m in model.metabolites:  
        cnt += 1
        charge = m.charge
        if charge==None or type(charge)!=int:
            if verbose: print(cnt, ':', m.id)
            found_mids.append(m.id)
            

    if found_mids == []:
        if verbose: print("No metabolite with missing charge attribute found.")
    return found_mids



def check_missing_formulas(model, verbose=True):
    """Check if all metabolites have a formula attribute.
    
    Args:
        model (cobra.Model): target model. 
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of matabolites missing the formula attribute.
    
    """ 
    
    found_mids = []
    cnt = 0
    for m in model.metabolites:  
        formula = m.formula
        if formula==None or formula=='':
            cnt += 1
            if verbose: print(cnt, ':', m.id)
            found_mids.append(m.id)
            

    if found_mids == []:
        if verbose: print("No metabolite with missing formula attribute found.")
    return found_mids



def check_artificial_atoms(model, preview=None, verbose=True):
    """Check if artificial atoms like 'R' and 'X' are present.
    
    Args:
        model (cobra.Model): target model. 
        preview (int): maximum number of metabolite IDs to show for each artificial atom. If ``None``, they will be all displayed.
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of matabolites with artificial atoms.
    """ 
    
    found_mids = []
    atom_to_mids = {}
    for m in model.metabolites:  
        formula = m.formula
        if formula == None or formula=='':
            continue  # there are dedicated functions for this
        
        
        # Matches any uppercase letter (A-Z) followed by zero or more lowercase letters (a-z)
        atoms = set(re.findall(r'[A-Z][a-z]*', formula))
        base = set(['C', 'H', 'O', 'N', 'P', 'S'])
        metals = set(['Ag', 'Fe', 'Co', 'As', 'Ca', 'Cd', 'Cl', 'Cu', 'Hg', 'K', 'Mg', 'Mo', 'Na', 'Ni', 'Se', 'Zn', 'Mn'])
        safe_atoms = base.union(metals)
        strange = atoms - safe_atoms
        
        
        # popupale te 'found_mids' list
        if len(strange) > 0: 
            found_mids.append(m.id)
        
        
        # iterate each 'strange' atom to create a dict: 
        for atom in strange: 
            if atom not in atom_to_mids.keys(): 
                atom_to_mids[atom] = []
            atom_to_mids[atom].append(f'{m.id}[{m.formula}]')
            
            
    # show the create 'atom_to_mids' dict: 
    for i, atom in enumerate(atom_to_mids.keys()): 
        mids_to_show = ', '.join(atom_to_mids[atom])
        if preview != None:
            if len(atom_to_mids[atom]) > preview:
                mids_to_show = f"{', '.join(atom_to_mids[atom][:preview])}, ... ({len(atom_to_mids[atom])} in total)"
        if verbose: print(i+1, ':', atom, ':', mids_to_show) 

        
    if found_mids == []:
        if verbose: print("No metabolite with artificial atoms found.")
    return found_mids



def get_unconstrained_bounds(model): 
    """Get the minimum and maximum bounds used in the model.
    
    Usually they correspond to the "unconstrained" negative and positive flux constants.
    
    Args:
        model (cobra.Model): target model. 
    
    Returns:
        tuple: min and max bounds.
    """
    un_lb, un_ub = 0, 0
    for r in model.reactions:
        if r.lower_bound < un_lb: un_lb = r.lower_bound
        if r.upper_bound > un_ub: un_ub = r.upper_bound
    return (un_lb, un_ub)



def reset_unconstrained_bounds(model): 
    """Set the uncontrained flux constants to 1000.
    
    Args:
        model (cobra.Model): target model. 
     
    """
    un_lb, un_ub = get_unconstrained_bounds(model)
    for r in model.reactions: 
        if r.lower_bound == un_lb: r.lower_bound = -1000
        if r.upper_bound == un_ub: r.upper_bound = 1000



def check_constrained_metabolic(model, verbose=True): 
    """Check the presence of constrained metabolic reactions.
    
    Metabolic reactions are here defined as those having more then 1 involved metabolites.
    Constrained reactions are here defined as those having bounds other then (0, 1000) or (-1000, 1000).
    
    Biomass assembly reactions are ignored.
    
    Args:
        model (cobra.Model): target model. 
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of reactions with constrained metabolic reactions.
    
    """
    
    found_rids = []
    cnt = 0
    biom_ids = search_biomass(model, verbose=False)
    for r in model.reactions: 
        if len(r.metabolites) == 1 or r.id in biom_ids: 
            continue # not interested in EXR/sink/demands, nor in biomass
            
        if r.bounds != (-1000, 1000) and r.bounds != (0, 1000): 
            cnt += 1
            if verbose: print(cnt, ':', r.id, ':', r.bounds, ':', r.reaction)
            found_rids.append(r.id)
    
    
    if found_rids == []:
        if verbose: print("No constrained metabolic reactions found.")
    return found_rids



def check_mass_unbalances(model, threshold=1e-5, verbose=True):
    """Check the presence of mass-unbalaned reactions in the model.
    
    EExchange, sink, demand, and biomass reactions are ignored.
    
    Args:
        model (cobra.Model): target model. 
        threshold (float): values below this treshold are considered as 0.
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of reactions with mass unbalanced.
    
    """
    
    found_rids = []
    cnt = 0
    biom_ids = search_biomass(model, verbose=False)
    for r in model.reactions:
        
        
        # exchage, sink, demand, and biomass reactions are excluded. 
        if len(r.metabolites)==1 or r.id in biom_ids: 
            continue 
        
        
        # get suggestions
        suggestions = r.check_mass_balance()
        
        # apply threshold: 
        suggestions = {atom: suggestions[atom] for atom in suggestions.keys() if suggestions[atom] >= threshold}
        
        # exclude charge is present:
        suggestions = {atom: suggestions[atom] for atom in list(set(list(suggestions.keys())) - set(['charge']))} 
        
        
        # log the suggestion:
        if len(suggestions.keys()) > 0:
            cnt += 1
            if verbose: print(cnt, ':', r.id, ':', r.reaction, ':', suggestions)
            found_rids.append(r.id)
    
    
    if found_rids == []:
        if verbose: print("No mass-unbalanced reactions found.")
    return found_rids



def check_charge_unbalances(model, threshold=1e-5, verbose=True):
    """Check the presence of charge-unbalaned reactions in the model.
    
    Exchange, sink, demand, and biomass reactions are ignored.
    
    Args:
        model (cobra.Model): target model. 
        threshold (float): values below this treshold are considered as 0.
        verbose (bool): if ``False``, don't print.
        
    Returns:
        list: IDs of reactions with charge unbalanced.
        
    """
    
    found_rids = []
    cnt = 0
    biom_ids = search_biomass(model, verbose=False)
    for r in model.reactions:
        
        
        # exchage, sink, demand, and biomass reactions are excluded. 
        if len(r.metabolites)==1 or r.id in biom_ids: 
            continue 
        
        
        # get suggestions
        suggestions = r.check_mass_balance()
        
        # retain only charge if present:
        suggestions = {atom: suggestions[atom] for atom in suggestions.keys() if atom == 'charge'} 
        
        # apply threshold: 
        suggestions = {atom: suggestions[atom] for atom in suggestions.keys() if suggestions[atom] >= threshold}
        
        
        # log the suggestion:
        if len(suggestions.keys()) > 0:
            cnt += 1
            if verbose: print(cnt, ':', r.id, ':', r.reaction, ':', suggestions)
            found_rids.append(r.id)
        
        
    if found_rids == []:
        if verbose: print("No charge-unbalanced reactions found.")
    return found_rids



def search_biomass(model, show_reaction=False, verbose=True):
    """Search for biomass reactions.
    
    Simple function involving just an exact match of substrings: 'biomass', 'growth', 'bof'.
    
    Args:
        model (cobra.Model): target model. 
        show_reaction (bool): whether to show also the reaction string
    
    """
    
    # define key objects: 
    substrings = ['biomass', 'growth', 'bof']
    found_rids = []
    cnt = 0
    
    
    def print_reaction(cnt, show_reaction, verbose):
        if verbose:
            if show_reaction: 
                print(cnt, ':', r.id, ':', r.name, ':', r.reaction, ':', r.bounds)
            else: 
                print(cnt, ':', r.id, ':', r.name, ':', r.bounds)
    
    
    # first try with reaction IDs:
    for r in model.reactions:
        for substring in substrings: 
            if substring.lower() in r.id.lower():
                cnt += 1
                print_reaction(cnt, show_reaction, verbose)
                found_rids.append(r.id)
    if found_rids != []:           
        return found_rids
                
        
    # then try with reaction names: 
    for r in model.reactions: 
        for substring in substrings: 
            if substring.lower() in r.name.lower():
                cnt += 1
                print_reaction(cnt, show_reaction, verbose)
                found_rids.append(r.id)
    if found_rids != []:           
        return found_rids
    
    
    # no biomass reactions found: 
    if found_rids == []:
        if verbose: print("No biomass reactions found.")
    return found_rids
    
    

def sanity_report(model):
    """Print a small sanity report.
    
    Get a sanity report calling some of the gempipe.curate functions.
    
    Args:
        model (cobra.Model): target model. 
    
    """
    
    version = importlib.metadata.metadata("gempipe")["Version"]
    print(f"gempipe v{version} - sanity_report")
    
    print("model ID:", model.id)
    
    G = len(model.genes)
    R = len(model.reactions)
    M = len(model.metabolites)
    uM = len(set([m.id.rsplit('_', 1)[0] for m in model.metabolites]))
    groups = len(model.groups)
    print('G:', G, 'R:', R, 'M:', M, 'uM:', uM, 'groups:', groups)
    
    comps = [c for c in model.compartments]
    print("Compartments:", sorted(comps))
    
    biom_ids = search_biomass(model, verbose=False)
    print("Biomass assemblies:", len(biom_ids), biom_ids)
    
    obj_ids = get_objectives(model)
    print("Objectives:", len(obj_ids), obj_ids)
    
    res = model.optimize()
    print("Optimization:", res.objective_value, res.status, f'({get_solver(model)})')
    
    print("Unconstrained LB-UB:", get_unconstrained_bounds(model))
    
    print("Bad EX_change notation:", len(check_exr_notation(model, verbose=False)))
    
    print("Sink/demand:", len(check_sink_demand(model, verbose=False)))
    
    print("Constrained metabolic:", len(check_constrained_metabolic(model, verbose=False)))
    
    print("With 'artificial' atoms:", len(check_artificial_atoms(model, verbose=False)))
    
    print(
        "Missing formulas - charges:",
        len(check_missing_formulas(model, verbose=False)), '-',
        len(check_missing_charges(model, verbose=False)),
    )
    
    print(
        "Mass - charge unbalances:",
        len(check_mass_unbalances(model, verbose=False)), '-',
        len(check_charge_unbalances(model, verbose=False)),
    )
    