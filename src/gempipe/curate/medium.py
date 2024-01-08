

def reset_uptakes(model):
    """Set all the EX_change reactions to (0, 1000).
    
    Args:
        model (cobra.Model): target model.

    """
    
    for r in model.reactions:
        if len(r.metabolites)==1 and r.id.startswith('EX_'):
            r.bounds = (0, 1000)
            
            
            
def set_bounded_uptakes(model, uptakes):
    """Set uptake of nutrients.
    
    Can be expressed as either concentrations [mmol/L] or fluxes [mmol/(h*gDW)], depending on your aims.
    
    Args:
        model (cobra.Model): target model.
        secretions (dict): dictionary keyed by EX_change reaction ID.

    """
    
    modeled_rids = [r.id for r in model.reactions]
    for ex_rid in uptakes.keys():
        if ex_rid not in modeled_rids: 
            raise Exception(f"{ex_rid} does not exists in this model.")
        model.reactions.get_by_id(ex_rid).lower_bound = -uptakes[ex_rid]
        
        
        
def set_bounded_secretions(model, secretions):
    """Set secretion of metabolites.
    
    Can be expressed as either concentrations [mmol/L] or fluxes [mmol/(h*gDW)], depending on your aims.
    
    Args:
        model (cobra.Model): target model.
        secretions (dict): dictionary keyed by EX_change reaction ID.

    """
    
    modeled_rids = [r.id for r in model.reactions]
    for ex_rid in secretions.keys():
        if ex_rid not in modeled_rids: 
            raise Exception(f"{ex_rid} does not exists in this model.")
        model.reactions.get_by_id(ex_rid).lower_bound = secretions[ex_rid]
        
        
        
def set_unbounded_exchanges(model, exrs):
    """Define which melecules require a free exchange (-1000, 1000).
    
    Args:
        model (cobra.Model): target model.
        exrs (list): IDs of EX_change reactions to unbound.

    """
    
    modeled_rids = [r.id for r in model.reactions]
    for ex_rid in exrs:
        if ex_rid not in modeled_rids: 
            raise Exception(f"{ex_rid} does not exists in this model.")
        model.reactions.get_by_id(ex_rid).bounds = (-1000, 1000)