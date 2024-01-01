import cobra



def solve_duplicate_m(logger, outdir):
    
    
    # log some message: 
    logger.info("Detecting duplicated metabolites using MetaNetX annotations...")
    
    
    # set up the cobra solver
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"
    
    
    # load the final draft panmodel
    draft_panmodel = cobra.io.load_json_model(outdir + 'draft_panmodel.json')
    
    
    logger.info(draft_panmodel.name)
    
    
    return 0