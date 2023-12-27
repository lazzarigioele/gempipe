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
        return f'({r.gene_reaction_rule}) or ({ref_r.gene_reaction_rule})'

    
    
def search_first_synonym(r, ref_model, remove_h=False ): 
    
    
    # check if the reference contains an equal reaction based on reactants and products mids. 
    for ref_r in ref_model.reactions:
        found_synonym = check_same_mids(r, ref_r, remove_h)
        if found_synonym: 
            return ref_r

        
    return None
    # pay attention if ATPM / NGAM is returned !
    
    
    
def get_detailed_reaction(r): 
    
    
    # make esplicit formula and charge after each metabolite
    reaction = r.reaction
    reaction = ' ' + reaction + ' '  # force each mid to be surrounded by spaces
    for m in r.metabolites:
        reaction = reaction.replace(f' {m.id} ', f' {m.id}[{m.formula},{m.charge}] ')
    reaction = reaction[1:-1]  # remove extreme spaces
    
    
    # include balancing information
    if r.check_mass_balance() != {}:
        reaction = '[[UNB]] ' + reaction
    
    return reaction 
    


def add_new_reaction(r, ref_model, addedms_logger):
    
    
    # Notes for manual correction of fatty-acid-ACP formulas: 
    # Take as example hexACP: 
    # C 390 H 613 O 143 N 96 P 1 S 3    (-1) <--- proposed
    # C 6 H 11 O X     (0) <--- desired 
    # C 384 H 602 O 142 N 96 P 1 S 3     <--- to subtract
    
    
    # get metabolites already modeled
    reference_mids = set([m.id for m in ref_model.metabolites])
    
    
    # first add missing metabolites:
    for m in r.metabolites: 
        if m.id not in reference_mids: 
            
            
            # create a new metabolite copying the attributes. 
            new_m = cobra.Metabolite(m.id)
            new_m.name = m.name
            
            # define formula and charge (here some manual corrections can be made)
            new_m.formula = m.formula
            new_m.charge = m.charge
            
            # actually add the metabolite
            ref_model.add_metabolites([new_m])
            new_m = ref_model.metabolites.get_by_id(m.id)
            
            # set compartment:
            new_m.compartment = new_m.id.rsplit('_', 1)[-1]   
            
            print(f'id="{new_m.id}"\tformula="{new_m.formula}"\tcharge="{new_m.charge}"\tname="{new_m.name}"', file=addedms_logger)
            
    
    # now add missing reaction: 
    # first create a new empty reaction
    new_r = cobra.Reaction(r.id)
    new_r.name = r.name
    
    # actually add the reaction
    ref_model.add_reactions([new_r])
    new_r = ref_model.reactions.get_by_id(r.id)
    
    # define the reaction (here some manual curations can be made)
    new_r_reaction = r.reaction
    new_r.build_reaction_from_string(new_r_reaction) 
    
    # set bounds
    new_r.bounds = r.bounds
    
    # copy GPR: 
    new_r.gene_reaction_rule = r.gene_reaction_rule 
    new_r.update_genes_from_gpr()
    
    

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
    addedms_logger = open('working/expansion/added_metabolites.txt', 'w')
    for r in draft_panmodel.reactions:
        if r.id == 'Growth': continue  # using reference biomass definition
        gpr = r.gene_reaction_rule
        reaction = get_detailed_reaction(r)
        
        
        # define key objects: 
        action = '-'
        same_rid = False
        same_mids = False
        same_fc = False
        same_gids = False
        ref_reaction = '-'
        final_reaction = '-'
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
            ref_reaction = get_detailed_reaction(ref_r)
            
            
            # check if reactants and products are the same: 
            same_mids = check_same_mids(r, ref_r, remove_h=True)
            if same_mids: 
                # check if also formula and charge correspond: 
                same_fc = check_same_fc(r, ref_model, remove_h=True)
                
                
            # check if the gpr is the same, and define final gpr: 
            same_gids = check_same_gids(r, ref_r)
            if not same_gids: 
                # if different, glue together the two gprs, and update the reference:  
                action = 'update_gpr'
                final_gpr = get_glued_gpr(r, ref_r)
                ref_r.gene_reaction_rule = gpr
                ref_r.update_genes_from_gpr()
            else: # if same gene set, just keep the definition from the reference 
                action = 'ignore'
                final_gpr = ref_gpr
                             
                
        else:  # new reaction to be included:
            action = 'add'
            add_new_reaction(r, ref_model, addedms_logger)
            new_r = ref_model.reactions.get_by_id(r.id)
            final_reaction = get_detailed_reaction(new_r)
         
        
        # populate the results dataframe: 
        new_row = {
            'rid': r.id, 'action': action, 
            'reaction': reaction, 'same_rid': same_rid, 'synonym': synonym, 'ref_reaction': ref_reaction,
            'same_mids': same_mids, 'same_fc': same_fc, 'final_reaction': final_reaction, 
            'same_gids': same_gids, 'gpr': gpr, 'ref_gpr': ref_gpr, 'final_gpr': final_gpr, }
        results_df.append(new_row)
    
    
    # save results dataframe to disk
    addedms_logger.close()
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/expansion/results.csv')
    
    
    # log some message
    draft_panmodel_exp = ref_model  # after the expansion
    logger.info(f"Done, {' '.join(['G:', str(len(draft_panmodel_exp.genes)), '|', 'R:', str(len(draft_panmodel_exp.reactions)), '|', 'M:', str(len(draft_panmodel_exp.metabolites))])}.")
    
    
    # count how many unbalanced reactions remained
    unbalanced = results_df[results_df['rid'].str.startswith('EX_')==False]
    unbalanced = results_df[(results_df['final_reaction'].str.startswith('[[UNB]] ')) | (results_df['ref_reaction'].str.startswith('[[UNB]] '))]
    logger.info(f"Terminated with {len(unbalanced)} unbalanced reactions.")
    
    
    # finally save the new draft panmodel
    cobra.io.save_json_model(draft_panmodel_exp, 'working/expansion/draft_panmodel.json')
    logger.debug("New draft pan-model saved to 'working/expansion/draft_panmodel.json'.")
    
    
    return 0

    #  G: 975 | R: 771 | M: 662.
    #  G: 1241 | R: 1447 | M: 1204. 258 unbalanced