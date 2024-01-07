

def verify_egc(model, mid, escher=False): 
    """Test the presence of energy-generating cycles (EGCs). 
    
    Can also output a model for Escher, with just the reactions composing the cycle. 
    
    Args:
        model: ``cobra.Model`` to verify. Must be encoded with the BiGG notation.
        mid: metabolite ID for which the EGC must be checked. Warning: must be without compartment, so for example ``atp`` instead of ``atp_c``. 
        escher: save a reduced ``cobra.Model`` in the current directory. To be loaded in Escher. 

    """
    
    # changes as not permament: 
    with model: 
        
        # close all the exchange reactions: 
        for r in model.reactions:
            if r.id.startswith('EX') and len(r.metabolites) ==1: 
                r.bounds = (0,1000)

                
        # create a dissipation reaction: 
        dissip = cobra.Reaction(f'__dissip__{mid}')
        model.add_reactions([dissip])
        dissip = model.reactions.get_by_id(f'__dissip__{mid}')
        
        
        # define the dissipation reaction:
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
            dissip_string = 'glu__L_c + h2o_c --> akg_c + nh4_c + 2.0 h_c'
        elif mid == 'q8h2':
            dissip_string = 'q8h2_c --> q8_c + 2.0 h_c'
        else: 
            raise Exception("Metabolite ID (mid) not recognized.") 
        dissip.build_reaction_from_string(dissip_string)
        
        
        # set the objective and optimize: 
        model.objective = f'__dissip__{mid}'
        res = model.optimize()
        
        
        # log some messages
        print("Testing the following dissipation reaction:", dissip.reaction)
        print("Objective value:", res.objective_value)
    
    
        # get suspect !=0 fluxes (if any)
        if res.objective_value != 0: 
            fluxes = res.fluxes
            fluxes_interesting = fluxes[fluxes != 0]
            print(fluxes_interesting)
            
            
            # create a model for escher
            if escher:  
                model_copy = model.copy()
                all_rids = [r.id for r in model_copy.reactions]
                to_delete = set(all_rids) - set(fluxes_interesting.index)
                model_copy.remove_reactions(to_delete)
                cobra.io.save_json_model(model_copy, f'__dissip__{mid}' + '.json')
                print(f'__dissip__{mid}' + '.json', "saved in current directory.")