import os


import pandas as pnd
import cobra



def get_first_occurrence(puremid, model):
    
    
    # get the fisrt metabolite in a model matching the provided 'puremid'.
    # 'puremid' is a metabolite ID without compartment. 
    for m in model.metabolites: 
        if puremid == m.id.rsplit('_', 1)[0] :
            return m
    return None



def translate_targets(logger, model, to_translate):
    # attemps to rewrite reactions in a model, in a way that duplicate metabolites are solved.

    
    # define key objects: 
    modeled_mids  = [m.id for m in model.metabolites]
    results_df = []  # tabular results
    
    
    # for each reaction in the model, solve each eventual duplicate
    for r in model.reactions:
        
        # get mids, puremids, compartements, for each involved metabolite. 
        r_mids     = [m.id                   for m in r.metabolites]
        r_puremids = [m.id.rsplit('_', 1)[0] for m in r.metabolites]
        r_comps  =   [m.id.rsplit('_', 1)[1] for m in r.metabolites]
        
        
        # define 'old' and future reaction strings
        reaction_old = r.reaction 
        reaction_new = r.reaction
        reaction_new = ' '+reaction_new+' '  # for each mid to be surrounded by spaces.
        
        
        # record new metabolites defined in order to translate this reaction
        added_mids = []
        # duplicate ~~> translated mids in this reaction
        translated_mids = []
        
        
        # iterate thrugh each duplicate-replacement pair: 
        for duplicate in to_translate.keys(): 
            replacement = to_translate[duplicate]
            

            # check if this reaction contains the metabolite to change: 
            if duplicate in r_puremids: 
                translated_mids.append(f'{duplicate} ~~> {replacement}')
                

                # iterate involved metabolites until we reach the one to change: 
                for mid, puremid, comp in zip(r_mids, r_puremids, r_comps): 
                    if puremid == duplicate: 
                        # check if the new metabolite already exists.
                        if f'{replacement}_{comp}' not in modeled_mids:

                            # retrive the metabolite to use as template (coming from another compartment)
                            good_m = get_first_occurrence(replacement, model)
                            # create new metabolite:
                            new_m = cobra.Metabolite(f'{replacement}_{comp}')
                            model.add_metabolites([new_m])
                            # copy formula, charge and annotation
                            new_m = model.metabolites.get_by_id(f'{replacement}_{comp}')
                            new_m.formula = good_m.formula
                            new_m.charge  = good_m.charge
                            new_m.annotation = good_m.annotation
                            # set the right compartment
                            new_m.compartment = comp

                            # update the 'modeled_mids'
                            modeled_mids.append(f'{replacement}_{comp}')
                            
                            # record new metabolites needed by this translation
                            added_mids.append(f'{replacement}_{comp}')
                            

                        # edit the reaction
                        reaction_new = reaction_new.replace(' '+mid+' ', ' '+f'{replacement}_{comp}'+' ')
                        
                        
        # remove the spaces
        reaction_new = reaction_new[1:-1]
        
        if reaction_new != reaction_old:  # if some changes were made: 
            # apply updated reaction
            r.build_reaction_from_string(reaction_new)


            # populate tabular results
            translated_mids =  f'{len(translated_mids)}: {translated_mids}'
            added_mids = f'{len(added_mids)}: {added_mids}'
            if added_mids == '0: []': added_mids = '-'
            bal_suggestions = str(r.check_mass_balance())
            if bal_suggestions == '{}': bal_suggestions = '-'
            if r.id.startswith('EX_') and len(r.metabolites)==1: bal_suggestions = '-'
            results_df.append({
                'rid': r.id, 'reaction_old': reaction_old, '~~>': '~~>', 'reaction_new': reaction_new, 
                'translated_mids': translated_mids, 'added_mids': added_mids, 'bal_suggestions': bal_suggestions})
    
    
    # save tabular results: 
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/duplicates/dup_m_translations.csv')



def sort_translation_dictionary(to_translate, model, refmodel, reffree):
    
    
    # 'to_translate' contains duplicated metabolites divided in groups. 
    # here, for each group, we select only one 'good' metabolite to maintain. 
    # Therefore, dictionary 1-to-many will be converted in a 1-to-1. 
    to_translate_11 = {}
    
    
    # get all mids without compartment
    puremids_refmodel = set([m.id.rsplit('_', 1)[0] for m in refmodel.metabolites])
    puremids_reffree = set([m.id.rsplit('_', 1)[0] for m in reffree.metabolites])
    
    
    results_df = []
    for group in to_translate.values():
        
        # determine the good metabolite to represent the group:
        good_puremid = list(group)[0]  # begin with the first
        for puremid in group:
            if puremid in puremids_refmodel: 
                good_puremid = puremid  # give precedence to refmodel
                break
                
                
        # solve in a 1-to-1 dict
        for puremid in group:
            if puremid != good_puremid: 
                # get formula, charge, and flags:
                
                
                # for this metabolite to replace
                m = get_first_occurrence(puremid, model)
                formula = m.formula
                charge = m.charge
                in_refmodel = 'RM' if puremid in puremids_refmodel else '--'
                in_reffree = 'RF' if puremid in puremids_reffree else '--'
                
                # for the replacement metabolite
                good_m = get_first_occurrence(good_puremid, model)
                good_formula = good_m.formula
                good_charge = good_m.charge
                good_in_refmodel = 'RM' if good_puremid in puremids_refmodel else '--'
                good_in_reffree = 'RF' if good_puremid in puremids_reffree else '--'
                
                
                # WARNING: by design, skip if both metabolite are ancoded in the refmodel:
                if in_refmodel == 'RM' and good_in_refmodel == 'RM':
                    continue
                
                
                # log the future edits: 
                new_row = {
                    'duplicated': puremid, 'd_formula': formula, 'd_charge': charge, 'd_in_refmodel': in_refmodel, 'd_in_reffree': in_reffree,
                    '~~>': '~~>',
                    'replacement': good_puremid, 'r_formula': good_formula, 'r_charge': good_charge, 'r_in_refmodel': good_in_refmodel, 'r_in_reffree': good_in_reffree}
                results_df.append(new_row)
                
                
                # populate the t-to-1 dictionary
                to_translate_11[puremid] = good_puremid
                
                
    # save tabular results
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/duplicates/dup_m_edits.csv')         
    return to_translate_11

    
    
def get_translation_dictionary_MNX(model): 

    
    # several notes on mnx: 
    # 1) some metabolites sometimes need to stay separated, eg nh3 and nh4. 
    # 2) even MNX is not perfect, for example lald__L and abt__L are clearly different species.
    # to_skip = ['nh3', 'nh4', 'lald__L' , 'abt__L']
    to_skip = ['lald__L' , 'abt__L', 'h2co3', 'hco3']
    
    
    # Dictionary grouping duplicated metabolites. The key is just one metabolite of the group. 
    # Later a 'good' member will be determined for each group. Here each group could have different f/c.
    to_translate = {}


    # get all mids without compartment
    puremids_model = set([m.id.rsplit('_', 1)[0] for m in model.metabolites])


    # put here metabolites already parsed (puremids):
    parsed = set()

    
    # iterate through all metabolites in this model
    for m in model.metabolites:

        # get the mid without compartement: 
        puremid = m.id.rsplit('_', 1)[0]
        if puremid in to_skip: continue

        # already parsed the same metabolite (probably from a different compartment):
        if puremid in parsed:
            continue
        parsed.add(puremid)

        # get the annotations
        try: biggannots = m.annotation['bigg.metabolite']
        except: continue

        
        # search if some of the other annotations are present in the model: 
        for biggid in biggannots: 
            if biggid != puremid and biggid in puremids_model:
        

                # already parsed the same metabolite (probably froma different compartment):
                if biggid in parsed:
                    continue
                parsed.add(biggid)

                
                # extract the corresponding metabolite in the model: 
                m2 = get_first_occurrence(biggid, model) 

                
                # populate the dictionary: 
                if puremid not in to_translate.keys():
                    to_translate[puremid] = set()
                    to_translate[puremid].add(puremid) # include the key in the group
                to_translate[puremid].add(m2.id.rsplit('_', 1)[0])
                    
                    
    return to_translate


def parse_disconnected_metabolites(logger, model):
    
    
    # remove metabolites involved to 0 reactions. 
    to_remove = []
    for m in model.metabolites:
        if len(m.reactions) == 0: 
            to_remove.append(m)
    model.remove_metabolites(to_remove)
    logger.debug(f"Removed {len(to_remove)} disconnected metabolites: {[m.id for m in to_remove]}.")

    
def check_duplicate_r_MNX(model, curated): 
    """
    Solve duplicate reactions based on their MNX annotation. 
    """
    
    # reactions (rids) contained in the curated model. 
    curated_rids = [r.id for r in curated.reactions]

    
    # PART A) build a dictionary keyed by MNX id, containing the 
    # corresponing set of duplicated reactions.
    mnx_dict = {}
    for r in model.reactions:

        if r.id.startswith("EX_") : 
            continue  # not interested in EX_ reaction

        # get the MNX annotations for this reaction: 
        try: mnx_annots = r.annotation['metanetx.reaction']
        except: mnx_annots = []
        if type(mnx_annots) == str: 
            mnx_annots = [mnx_annots]
        for mnx_id in mnx_annots: 
            if mnx_id not in mnx_dict.keys():
                mnx_dict[mnx_id] = set() 
            mnx_dict[mnx_id].add(r.id)
            
    
    # PART B) iterate the dictionary of duplicates and solve duplications.
    to_delete = [] # reactions to delete
    for mnx_id in mnx_dict.keys():
        if len(mnx_dict[mnx_id] ) ==1: 
            continue # interested only in duplicated reactions. 
        # print the rids of these 'repeating group'
        print(mnx_dict[mnx_id])


        # STEP 1) print their reactions, and identify the 'cu' (curated) member (if any). 
        cu_rid = None
        for rid in mnx_dict[mnx_id]:
            r = model.reactions.get_by_id(rid)
            msg = 'cu' if rid in curated_rids else '--'
            if msg=='cu': cu_rid = rid
            print('.     ', msg, rid, r.reaction, r.bounds)


        # STEP 2) check if they are really the same in terms of reactants / products.
        # Take the first item as reference: 
        rid0 = list(mnx_dict[mnx_id])[0]
        r0 = model.reactions.get_by_id(rid0)
        reacs0 = [m.id for m in r0.reactants]
        prods0 = [m.id for m in r0.products]
        # Iterate the others:
        same_ids = True
        for rid in list(mnx_dict[mnx_id])[1:]:
            r = model.reactions.get_by_id(rid)
            reacs = [m.id for m in r.reactants]
            prods = [m.id for m in r.products]
            if set(reacs0) == set(reacs) and set(prods0) == set(prods):
                continue
            elif set(reacs0) == set(prods) and set(prods0) == set(reacs):
                continue
            else:
                same_ids = False
                break
        print('_______Same metabolites?', same_ids)


        # STEP 3) check if they have the same GPR
        gpr0 = r0.gene_reaction_rule
        same_gpr= True
        print('.     ', gpr0)
        for rid in list(mnx_dict[mnx_id])[1:]:
            r = model.reactions.get_by_id(rid)
            gpr = r.gene_reaction_rule
            print('.     ', gpr)
            if gpr == gpr0: 
                continue
            else:
                same_gpr = False
        print('_______Same GPR?', same_gpr)
        
        
        # STEP 4) Keeping just 1 of the group.
        # The GPR of the 'keeping_rid' will be updated, concatenating using ORs.
        if cu_rid != None: keeping_rid = cu_rid
        else: 
            # take the first reaction not having all metabolites in _e compartment.
            # (assuming that at least 1 'good' reaction exists. 
            for rid in mnx_dict[mnx_id]:
                r = model.reactions.get_by_id(rid)
                if not all([True if m.id.endswith('_e') else False for m in r.metabolites ]):
                    keeping_rid = rid
                    break
        keeping_r = model.reactions.get_by_id(keeping_rid)
        keeping_gpr = keeping_r.gene_reaction_rule
        for rid in mnx_dict[mnx_id]: 
            if rid != keeping_rid:
                r = model.reactions.get_by_id(rid)
                to_delete.append(r)
                if keeping_gpr == '': keeping_gpr = r.gene_reaction_rule
                elif r.gene_reaction_rule == '': pass
                else: keeping_gpr = keeping_gpr + ' or ' + r.gene_reaction_rule
                print('======>', rid, 'marked for DELETION!')
        keeping_r.gene_reaction_rule = keeping_gpr
        keeping_r.update_genes_from_gpr()
                
        
        print() 
    model.remove_reactions(to_delete)
            
    

def solve_duplicates(logger, outdir,  identity, coverage, refmodel):
    
    
    # log some message: 
    logger.info("Detecting duplicated metabolites using MetaNetX annotations...")
    
    
    # create subdirs without overwriting: 
    os.makedirs('working/duplicates/', exist_ok=True)
    
    
    # set up the cobra solver
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"
    
    
    # load the needed models: 
    draft_panmodel = cobra.io.load_json_model(outdir + 'draft_panmodel.json')
    reffree = cobra.io.load_json_model(f'working/free/draft_panmodel_{identity}_{coverage}.json')
    refmodel_basename = os.path.basename(refmodel)
    refmodel = cobra.io.load_json_model(f'working/brh/{refmodel_basename}.refmodel_translated.json')
    
    
    # get the translation dictionary:                                       
    to_translate = get_translation_dictionary_MNX(draft_panmodel)  # 1-to-many
    to_translate = sort_translation_dictionary(to_translate, draft_panmodel, refmodel, reffree)  # 1-to-1
                      
                      
    # translate reactions with the 'good' metabolites; 
    translate_targets(logger, draft_panmodel, to_translate)
    
    
    # remove disconnected metabolites:
    parse_disconnected_metabolites(logger, draft_panmodel)
    
    
    # remove duplicate reaction
    check_duplicate_r_MNX(draft_panmodel, refmodel)
    
    
    # save deduplicated model
    cobra.io.save_json_model(draft_panmodel, f'working/duplicates/draft_panmodel_{identity}_{coverage}_ddM.json')
    
                             
    return 0