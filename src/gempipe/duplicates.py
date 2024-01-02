import os


import pandas as pnd
import cobra



def get_first_occurrence(puremid, model):
    for m in model.metabolites: 
        if puremid == m.id.rsplit('_', 1)[0] :
            return m
    return None



def translate_targets(model, to_translate):
    # attemps to rewrite reactions in a model, in a way that duplicate metabolites are solved.

    
    modeled_mids  = [m.id for m in model.metabolites]
    r_to_delete = []
    
    
    results_df = []
    for duplicate in to_translate.keys(): 
        replacement = to_translate[duplicate]
            
        for r in model.reactions: 

            r_mids     = [m.id                   for m in r.metabolites]
            r_puremids = [m.id.rsplit('_', 1)[0] for m in r.metabolites]
            r_comps  =   [m.id.rsplit('_', 1)[1] for m in r.metabolites]

            # check if this reaction contains the metabolite to change: 
            if duplicate in r_puremids: 

                reaction_old = r.reaction

                # put spaces at the beginning and end of the reaction, in order to
                # later use the str.replace() built-in function.  
                reaction_new = r.reaction
                reaction_new = ' ' + reaction_new + ' '

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
                            added = model.metabolites.get_by_id(f'{replacement}_{comp}')
                            added.formula = good_m.formula
                            added.charge  = good_m.charge
                            added.annotation = good_m.annotation
                            # set the right compartment
                            added.compartment = comp

                            # update the list
                            modeled_mids.append(f'{replacement}_{comp}')

                            #print(f'ADDED {replacement}_{comp} FOR {duplicate}')

                        # edit the reaction
                        reaction_new = reaction_new.replace(' '+mid+' ', ' '+f'{replacement}_{comp}'+' ')
                # remove the spaces
                reaction_new = reaction_new[1:-1]
                # apply updated reaction
                r.build_reaction_from_string(reaction_new)
                
                results_df.append({'duplicate': duplicate, 'replacement': replacement, 'reaction_old': reaction_old, '~~>': '~~>', 'reaction_new': reaction_new})

                #print(f'EDITED {r.id} AS {r.reaction} FOR {duplicate}')

                # in some cases, an empty reaction could be created:
                if r.reaction == ' <=> ' or r.reaction == ' --> ' or r.reaction == ' <-- ':
                    #print('.   ', 'marked for deletion!')
                    r_to_delete.append(r)

                if r.check_mass_balance() != {} and len(r.metabolites) != 1: 
                    #print("UNBALANCED", r.check_mass_balance())
                    pass
    
    # remove flagged reactions
    model.remove_reactions(r_to_delete)
    
    # save tabular results: 
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/duplicates/dup_m_translations.csv')



def sort_translation_dictionary(to_translate, model, refmodel, reffree):
    # 'to_translate' contains duplicated metabolites divided in groups. 
    # here, for each group, we keep one 'good' metabolite to maintain. 
    # dictionary 1-to-many will be converted in a 1-to-1. 
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
                good_puremid = puremid  # give precedente to refmodel
                break
                
        # solve in a 1-to-1 dict
        for puremid in group:
            if puremid != good_puremid: 
                   
                # get formula and charge to log
                # for this metabolite to replace: 
                m = get_first_occurrence(puremid, model)
                formula = m.formula
                charge = m.charge
                # for the replacement metabolite: 
                good_m = get_first_occurrence(good_puremid, model)
                good_formula = good_m.formula
                good_charge = good_m.charge
                
                # determine some flags for the presence in models
                # for this metabolite to replace:
                in_refmodel = 'RM' if puremid in puremids_refmodel else '--'
                in_reffree = 'RF' if puremid in puremids_reffree else '--'
                # for the replacement metabolite: 
                good_in_refmodel = 'RM' if good_puremid in puremids_refmodel else '--'
                good_in_reffree = 'RF' if good_puremid in puremids_reffree else '--'
                
                
                # WARNING: by design, skip if both metabolite are ancoded in the refmodel
                if in_refmodel == 'RM' and good_in_refmodel == 'RM':
                    continue
                
                # log the future edits: 
                new_row = {
                    'duplicated': puremid, 'd_formula': formula, 'd_charge': charge, 'd_in_refmodel': in_refmodel, 'd_in_reffree': in_reffree,
                    '~~>': '~~>',
                    'replacement': good_puremid, 'r_formula': good_formula, 'r_charge': good_charge, 'r_in_refmodel': good_in_refmodel, 'r_in_reffree': good_in_reffree}
                results_df.append(new_row)
                
                # populate the dictionary
                to_translate_11[puremid] = good_puremid
                
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/duplicates/dup_m_edits.csv')         
    return to_translate_11

    
    
def get_duplicated_MNX(model): 

    
    # # Some notes: 
    # 1) some metabolites need to stay separated, eg nh3 and nh4. 
    # 2) even MNX is not perfect, for example lald__L and abt__L are clearly different species.
    # to_skip = ['nh3', 'nh4', 'lald__L' , 'abt__L']
    to_skip = []
    
    
    # Dictionary grouping duplicated metabolites. The key is just one metabolite of the group. 
    # Later a 'good' member will be determined for each group. Here each group could have different f/c.
    to_translate = {}


    # get all mids without compartment
    puremids_model = set([m.id.rsplit('_', 1)[0] for m in model.metabolites])


    # put here metabolites already parsed (pure mids):
    parsed = set()

    
    # iterate through all metabolites in this model
    cnt = 0
    for m in model.metabolites:

        # get the mid without compartement: 
        puremid = m.id.rsplit('_', 1)[0]

        # already parsed the same metabolite (probably from a different compartment):
        if puremid in parsed:
            continue
        parsed.add(puremid)

        # get the annotations
        try: biggannots = m.annotation['bigg.metabolite']
        except: continue

        
        # search if some of the annotations are present in the model: 
        for biggid in biggannots:
            # if this model contains other equivalent metabolites: 
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



def solve_duplicate_m(logger, outdir,  identity, coverage, refmodel):
    
    
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
    to_translate = get_duplicated_MNX(draft_panmodel)  # 1-to-many
    to_translate = sort_translation_dictionary(to_translate, draft_panmodel, refmodel, reffree)  # 1-to-1
                      
                      
    # translate reactions with the 'good' metabolites; 
    translate_targets(draft_panmodel, to_translate)
    
    
    return 0