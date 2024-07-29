# experimental code, not used in production.


from optlang.glpk_interface import Model as Problem
import cobra



def test_manual_FBA(model, obid='biomass_LPL60'):
    
    
    # get the stoichiometric matrix (S): rows:metabolites, columns:reactions, cells:stoichiometric coefficients.
    # metabolite and reaction IDs are assumed unique (that is, no reapeated row ID or col ID).
    S = cobra.util.array.create_stoichiometric_matrix(model, array_type='DataFrame')
    
    
    # define the problem description in json format (dict)
    problem_json = {}  
    problem_json['name'] = 'manual_FBA'
    problem_json['variables'] = []
    problem_json['constraints'] = []
    problem_json['objective'] = {}
    
    
    # some convertions to speed up the table accession
    rid_to_index = {rid: index for index, rid in enumerate(S.columns)}
    mid_to_index = {mid: index for index, mid in enumerate(S.index)}
    S = S.values  # convert to numpy.ndarray  (faster cell accession)


    # define Variables, that are fluxes through reactions:
    # (variables are constrained with a LB and an UB)
    for rid in rid_to_index.keys():
        r = model.reactions.get_by_id(rid)
        problem_json['variables'].append({'name': rid, 'lb': r.lower_bound, 'ub': r.upper_bound, 'type': 'continuous'})
        
    
    # define Constraints, that are mass balances for each metabolite:
    # (in steady state, each mass balance equals to 0)
    for mid, index_r in mid_to_index.items():  # iterate metabolites
        add_args = []
        
        
        for rid, index_c in rid_to_index.items():  # iterate reactions
            stoich_coeff = S[index_r, index_c]   # get the stoichiometric coefficient from S.
            if stoich_coeff == 0: continue   # would be a useless term (save time).
            
            
            # define the addend (contribution of this reaction to the current mid):
            add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': stoich_coeff}, {'type': 'Symbol', 'name': rid}]})
        # add final mass balance for this mid:
        problem_json['constraints'].append({'name': f'{mid}_constraint', 'expression': {'type': 'Add', 'args': add_args}, 'lb': 0, 'ub': 0, 'indicator_variable': None, 'active_when': 1})
    
        
    # finally define the objective (1*obid): 
    obdir = 'max'  # minimize or maximize
    problem_json['objective'] = {'name': 'my_objective', 'expression': {'type': 'Mul', 'args': [{'type': 'Number', 'value': 1.0}, {'type': 'Symbol', 'name': obid}]}, 'direction': obdir}
    
    
    # convert json/dict to optlang problem:
    problem = Problem.from_json(problem_json)
    status = problem.optimize()
    return problem.objective.value



def test_prio_gapfilling(model, Rset_dataframe, biom_id='Growth', minflux=0.01, ):
    
    
    # get the stoichiometric matrix (S): rows:metabolites, columns:reactions, cells:stoichiometric coefficients.
    # metabolite and reaction IDs are assumed unique (that is, no reapeated row ID or col ID).
    S = cobra.util.array.create_stoichiometric_matrix(model, array_type='DataFrame')
    
    
    # define the problem description in json format (dict)
    problem_json = {}  
    problem_json['name'] = 'prio_gapfilling'
    problem_json['variables'] = []
    problem_json['constraints'] = []
    problem_json['objective'] = {}
    
    
    # some convertions to speed up the table accession
    rid_to_index = {rid: index for index, rid in enumerate(S.columns)}
    mid_to_index = {mid: index for index, mid in enumerate(S.index)}
    S = S.values  # convert to numpy.ndarray  (faster cell accession)
    
    
    # define the metabolic reactions that are not part of the draft panmodel:
    excluding = set(Rset_dataframe.index.to_list())


    # define Variables, that are fluxes through reactions:
    # (variables are constrained with a LB and an UB)
    for rid in rid_to_index.keys():
        r = model.reactions.get_by_id(rid)
        lb = r.lower_bound
        ub = r.upper_bound
        if rid == biom_id: lb = minflux  # apply minflux for biomass formation
        problem_json['variables'].append({'name': rid, 'lb': lb, 'ub': ub, 'type': 'continuous'})
        if rid in excluding:  # in addition to the flux variable, add another binar variable (indicator variable).
            problem_json['variables'].append({'name': f'binary_{rid}', 'lb': 0, 'ub': 1, 'type': 'binary'})
        

    # define Constraints, that are mass balances for each metabolite:
    # (in steady state, each mass balance equals to 0)
    for mid, index_r in mid_to_index.items():  # iterate metabolites
        add_args = []
        
        
        for rid, index_c in rid_to_index.items():  # iterate reactions
            stoich_coeff = S[index_r, index_c]   # get the stoichiometric coefficient from S.
            if stoich_coeff == 0: continue   # would be a useless term (save time).
            
            
            # define the addend (contribution of this reaction to the current mid):
            add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': stoich_coeff}, {'type': 'Symbol', 'name': rid}]})
        # add final mass balance for this mid:
        problem_json['constraints'].append({'name': f'{mid}_constraint', 'expression': {'type': 'Add', 'args': add_args}, 'lb': 0, 'ub': 0, 'indicator_variable': None, 'active_when': 1})
   
        
    # Define additional Constraints, using the indicator variables.
    # Contraints defined here: yi * LBi <= vi <= yi * UBi  (for each reaction i in 'excluding')
    # Decompose the two terms: 
    # yi * LBi <= vi   ===>   vi - yi * LBi >= 0    
    # vi <= yi * UBi   ===>   vi - yi * UBi <= 0
    # Where: --- yi is the indicator variable for the reaction i
    #        --- LBi, UBi are lower- and upper bounds for the reaction i
    #        --- vi is the flux variable through reaction i 
    #        --- i are potential gap-filler, ie reactions that are not inlcuded in the draft.
    for rid, row in Rset_dataframe.iterrows(): 
        
        
        # vi - yi * LBi >= 0    ===dissecting===>  [vi]  [- yi * LBi]  [>= 0]
        add_args = []
        add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': 1}, {'type': 'Symbol', 'name': rid}]})  #  vi
        add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': -row['lb']}, {'type': 'Symbol', 'name': f'binary_{rid}'}]})  # - yi * LBi
        problem_json['constraints'].append({'name': f'indicator_{rid}_lb', 'expression': {'type': 'Add', 'args': add_args}, 'lb': 0, 'ub': 1000, 'indicator_variable': None, 'active_when': 1})  # >= 0
        
        
        # vi - yi * UBi <= 0    ===dissecting===>  [vi]  [- yi * UBi]  [<= 0]
        add_args = []
        add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': 1}, {'type': 'Symbol', 'name': rid}]})  #  vi
        add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': -row['ub']}, {'type': 'Symbol', 'name': f'binary_{rid}'}]})  # - yi * UBi
        problem_json['constraints'].append({'name': f'indicator_{rid}_ub', 'expression': {'type': 'Add', 'args': add_args}, 'lb': -1000, 'ub': 0, 'indicator_variable': None, 'active_when': 1})  # <= 0
    
        
    # finally define the objective :  sum_i ( (1/(1+n_score))*yi )
    obdir = 'min'  # minimize or maximize
    add_args = []
    for rid, row in Rset_dataframe.iterrows(): 
        penalty = 1 / (1 + row['n_score'])  # inversely proportional to the normalized score
        add_args.append({'type': 'Mul', 'args': [{'type': 'Number', 'value': penalty}, {'type': 'Symbol', 'name': f'binary_{rid}'}]}) 
        
        
    problem_json['objective'] = {'name': 'my_objective', 'expression': {'type': 'Add', 'args': add_args}, 'direction': obdir}
    
    
    # convert json/dict to optlang problem:
    problem = Problem.from_json(problem_json)
    return problem

    """
    status = problem.optimize()
    return problem.objective.value
    """