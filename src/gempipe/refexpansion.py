import os


import pandas as pnd
import cobra



def check_same_mids(r, ref_r, remove_h=False): 
    
    
    # get metabolites and products of each reaction: 
    reacs = [m.id for m in r.reactants]
    prods = [m.id for m in r.products]
    ref_reacs = [m.id for m in ref_r.reactants]
    ref_prods = [m.id for m in ref_r.products]
    
    
    # if requested, do not consider protons:
    if remove_h: 
        reacs  = [mid for mid in reacs if mid.rsplit('_', 1)[0] != 'h']
        prods  = [mid for mid in prods if mid.rsplit('_', 1)[0] != 'h']
        ref_reacs  = [mid for mid in ref_reacs if mid.rsplit('_', 1)[0] != 'h']
        ref_prods  = [mid for mid in ref_prods if mid.rsplit('_', 1)[0] != 'h']
        

    # test both ways as reactions could be written upside down
    if   set(reacs)==set(ref_reacs) and set(prods)==set(ref_prods): 
        return True
    elif set(reacs)==set(ref_prods) and set(prods)==set(ref_reacs): 
        return True
    else: 
        return False
    
    
    
def check_same_fc(r, ref_model, remove_h=False):

    
    # parse each metabolite of this reaction
    comparison = []
    for m in r.metabolites: 
        
        
        # if requested, do not consider protons:
        if remove_h:
            if m.id.rsplit('_', 1)[0] == 'h':
                continue
        
        
        # check if same molecular formula and charge: 
        ref_m = ref_model.metabolites.get_by_id(m.id)
        if (m.formula == ref_m.formula) and (m.charge == ref_m.charge):
            comparison.append(True)
        else: 
            comparison.append(False)
    return all(comparison)



def check_same_gids(r, ref_r):

    
    # get the gene sets
    gids = set([g.id for g in r.genes])
    ref_gids = set([g.id for g in ref_r.genes])
    
    
    # check if two reactions have the same gene set in their gprs.
    if gids == ref_gids: 
        return True
    else: return False



def get_glued_gpr(r, ref_r):
    
    
    # get the gene sets
    gids = set([g.id for g in r.genes])
    ref_gids = set([g.id for g in ref_r.genes])
    
    
    # glue together the gprs from two different reactions: 
    if len(gids) == 0 and len(ref_gids) == 0:
        return ''
    elif len(gids) == 0: 
        return ref_r.gene_reaction_rule
    elif len(ref_gids) == 0:
        return r.gene_reaction_rule
    else: 
        return f'({ref_r.gene_reaction_rule}) or ({r.gene_reaction_rule})'

    
    
def search_first_synonym(r, ref_model, remove_h=False ): 
    
    
    # check if the reference contains an equal reaction based on reactants and products mids. 
    for ref_r in ref_model.reactions:
        found_synonym = check_same_mids(r, ref_r, remove_h)
        if found_synonym: 
            return ref_r

        
    return None
    # pay attention if ATPM / NGAM is returned !
    


def ref_expansion(logger, identity, coverage): 
    
    
    
    # set up the cobra solver
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/expansion/', exist_ok=True)
    
    
    # log some messages
    logger.info("Expanding the reference model with new reactions taken from the reference-free reconstruction...")
    
    
    # load the reference and the reference-free reconstruction.
    ref_model = cobra.io.load_json_model('working/brh/ref_model_t.json')
    draft_panmodel = cobra.io.load_json_model(f'working/panmodel/draft_panmodel_{identity}_{coverage}.json')
    
    
    # begin the addition of new reactions to the ref_model, to form the final draft pan-model.
    reference_rids = [r.id for r in ref_model.reactions]
    results_df = []  # summarizing all the additions
    for r in draft_panmodel.reactions:
        gpr = r.gene_reaction_rule
        
        
        # define key objects: 
        same_rid = False
        same_mids = False
        same_fc = False
        same_gids = False
        ref_gpr = '-'
        final_gpr = '-'
        synonym = '-'
        
        
        # check if this rid is already modeled: 
        same_rid = r.id in reference_rids
        if same_rid:
            # get the reference reaction and gpr:
            ref_r = ref_model.reactions.get_by_id(r.id)
            
            
        else: # this rid is not modeled yet
            # search for synonyms in ref_model and take the first:
            ref_r = search_first_synonym(r, ref_model, remove_h=True)
            if ref_r != None: 
                synonym = ref_r.id
        
        
        # manage the reference reaction if any:
        if ref_r != None:
            ref_gpr = ref_r.gene_reaction_rule
            
            
            # check if reactants and products are the same: 
            same_mids = check_same_mids(r, ref_r, remove_h=True)
            if same_mids: 
                # check if also formula and charge correspond: 
                same_fc = check_same_fc(r, ref_model, remove_h=True)
                
                
            # check if the gpr is the same, and define final gpr: 
            same_gids = check_same_gids(r, ref_r)
            if not same_gids: 
                # if different, glue together the two gprs, and update the reference:  
                final_gpr = get_glued_gpr(r, ref_r)
                ref_r.gene_reaction_rule = gpr
                ref_r.update_genes_from_gpr()
            else: # if same gene set, just keep the definition from the reference 
                final_gpr = ref_gpr
                
                
        
            
            
            
            
        new_row = {
            'rid': r.id, 'same_rid': same_rid, 'synonym': synonym,
            'same_mids': same_mids, 'same_fc': same_fc, 
            'same_gids': same_gids, 'gpr': gpr, 'ref_gpr': ref_gpr, 'final_gpr': final_gpr}
        results_df.append(new_row)
    
    
    # save results dataframe to disk
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/expansion/results.csv')
    
    
    # log some message
    draft_panmodel_exp = ref_model  # after the expansion
    logger.info(f"Done, {' '.join(['G:', str(len(draft_panmodel_exp.genes)), '|', 'R:', str(len(draft_panmodel_exp.reactions)), '|', 'M:', str(len(draft_panmodel_exp.metabolites))])}.")
    
    
    # finally save the new draft panmodel
    cobra.io.save_json_model(draft_panmodel_exp, 'working/expansion/draft_panmodel.json')
    logger.debug("New draft pan-model saved to 'working/expansion/draft_panmodel.json'.")
    
    
    return 0