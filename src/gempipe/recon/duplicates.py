import os
import pickle


import pandas as pnd
import cobra


from ..commons import get_md5_string
from .refexpansion import mancor_to_dict



def get_first_occurrence(puremid, model):
    
    
    # get the fisrt metabolite in a model matching the provided 'puremid'.
    # 'puremid' is a metabolite ID without compartment. 
    for m in model.metabolites: 
        if puremid == m.id.rsplit('_', 1)[0] :
            return m
    return None


    
def get_translation_dictionary_mnx(model, mrmode='m'): 

    
    # even MNX is not perfect, for example lald__L and abt__L are clearly different species.
    to_skip_m = ['MNXM732880']  # 'MNXM732880': {'lald__L', 'abt__L'}
    to_skip_r = []
    if   mrmode == 'm': to_skip = to_skip_m
    elif mrmode == 'r': to_skip = to_skip_r
    
    
    # Dictionary grouping duplicated m/r. Later a 'good' member will be determined for each group. 
    # Here each group of metabolites could have different f/c.
    mnx_to_mrids = {}  # future 'to_translate'
    
    
    # iterate through all m/r in this model
    if   mrmode == 'm': items = model.metabolites
    elif mrmode == 'r': items = model.reactions
    for item in items:

        
        # get the annotations
        try: # mnxannots can be list of str (if just 1 element)
            if   mrmode == 'm': mnxannots = item.annotation['metanetx.chemical']
            elif mrmode == 'r': mnxannots = item.annotation['metanetx.reaction']
            if type(mnxannots) == str: mnxannots = [mnxannots]
        except: continue
        
        
        # populate the dict
        for mnxannot in mnxannots: 
            if mnxannot in to_skip: 
                continue
            if mnxannot not in mnx_to_mrids.keys():
                mnx_to_mrids[mnxannot] = set()
            if   mrmode == 'm': mnx_to_mrids[mnxannot].add(item.id.rsplit('_', 1)[0])  # without compartment. 
            elif mrmode == 'r': 
                if item.id.startswith("EX_") and len(item.metabolites)==1: 
                    continue  # by design not interest in EX_ reactions
                mnx_to_mrids[mnxannot].add(item.id) 
            
            
    # filter to keep only groups with >1 elements.
    to_translate = {}
    for mnxannot in mnx_to_mrids.keys(): 
        if len(mnx_to_mrids[mnxannot]) >= 2:
            to_translate[mnxannot] = mnx_to_mrids[mnxannot]
                    
                    
    return to_translate

    

def sort_translation_dictionary(to_translate, model, refmodel, reffree, mrmode='m'):
    
    
    # 'to_translate' contains duplicated metabolites divided in groups. 
    # here, for each group, we select only one 'good' metabolite to maintain. 
    # Therefore, dictionary 1-to-many will be converted in a 1-to-1. 
    to_translate_11 = {}
    
    
    # get all mids without compartment
    if   mrmode == 'm': 
        mrids_refmodel = set([m.id.rsplit('_', 1)[0] for m in refmodel.metabolites])
        mrids_reffree = set([m.id.rsplit('_', 1)[0] for m in reffree.metabolites])
    elif mrmode == 'r': 
        mrids_refmodel = set([r.id for r in refmodel.reactions])
        mrids_reffree = set([r.id for r in reffree.reactions])
    
    
    results_df = []
    for group in to_translate.values():
        
        # determine the good metabolite to represent the group:
        good_mrid = list(group)[0]  # begin with the first
        for mrid in group:
            if mrid in mrids_refmodel: 
                good_mrid = mrid  # give precedence to refmodel
                break
                
                
        # solve in a 1-to-1 dict
        for mrid in group:
            if mrid != good_mrid: 
                # get formula, charge, and flags:
                
                
                # for metabolties:  
                if mrmode == 'm':
                    
                    # for this metabolite to replace
                    m = get_first_occurrence(mrid, model)
                    formula = m.formula
                    charge = m.charge
                    in_refmodel = mrid in mrids_refmodel
                    in_reffree = mrid in mrids_reffree

                    # for the replacement metabolite
                    good_m = get_first_occurrence(good_mrid, model)
                    good_formula = good_m.formula
                    good_charge = good_m.charge
                    good_in_refmodel = good_mrid in mrids_refmodel
                    good_in_reffree = good_mrid in mrids_reffree


                    # WARNING: by design, skip if both metabolite are ancoded in the refmodel:
                    if in_refmodel and good_in_refmodel: 
                        continue
                        

                    # log the future edits: 
                    new_row = {
                        'duplicated': mrid, 'd_formula': formula, 'd_charge': charge, 'd_in_refmodel': in_refmodel, 'd_in_reffree': in_reffree,
                        '~~>': '~~>',
                        'replacement': good_mrid, 'r_formula': good_formula, 'r_charge': good_charge, 'r_in_refmodel': good_in_refmodel, 'r_in_reffree': good_in_reffree}
                    results_df.append(new_row)


                    # populate the t-to-1 dictionary
                    to_translate_11[mrid] = good_mrid
                
                
                # for reactions
                elif mrmode == 'r': 
                    
                    # for this reaction to replace
                    r = model.reactions.get_by_id(mrid)
                    in_refmodel = mrid in mrids_refmodel
                    in_reffree = mrid in mrids_reffree
                    mids_involved = set([m.id for m in r.metabolites if m.id.rsplit('_', 1)[0] != 'h'])  # exclude protons.
                    gids_involved = set([g.id for g in r.genes])  
                    
                    # for the replacement reaction
                    good_r = model.reactions.get_by_id(good_mrid)
                    good_in_refmodel = good_mrid in mrids_refmodel
                    good_in_reffree = good_mrid in mrids_reffree
                    good_mids_involved = set([m.id for m in good_r.metabolites if m.id.rsplit('_', 1)[0] != 'h'])  # exclude protons.
                    good_gids_involved = set([g.id for g in good_r.genes])  
                    
                    
                    # determine if exaclty same reactants (except for protons).
                    # WARNING: this way we exclude transporters involved in different compartments
                    if mids_involved != good_mids_involved:
                        continue
                    
                    
                    # WARNING: by design, skip if both reactions are ancoded in the refmodel:
                    if in_refmodel and good_in_refmodel:
                        continue
                        
                        
                    # WARNING: by design, if two equivalment reactions have different gene sets, glue the gprs:
                    same_geneset = gids_involved == good_gids_involved
                        
                        
                    # log the future edits: 
                    new_row = {
                        'duplicated': mrid, 'd_reaction': r.reaction, 'd_gpr': r.gene_reaction_rule, 'd_in_refmodel': in_refmodel, 'd_in_reffree': in_reffree,
                        '~~>': '~~>',
                        'replacement': good_mrid, 'r_reaction': good_r.reaction, 'r_gpr': good_r.gene_reaction_rule, 'r_in_refmodel': good_in_refmodel, 'r_in_reffree': good_in_reffree, 'same_geneset': same_geneset}
                    results_df.append(new_row)


                    # populate the t-to-1 dictionary
                    to_translate_11[mrid] = good_mrid
                
     
    # produce the final results table: 
    results_df = pnd.DataFrame.from_records(results_df)
    
    
    # erease meaningless fileds when not using a reference
    if refmodel.id == '__EMPTY__': 
        results_df['d_in_refmodel'] = '-'
        results_df['r_in_refmodel'] = '-'
        
        
    # save tabular results
    if   mrmode=='m': results_df.to_csv('working/duplicates/dup_m_edits.csv')  
    elif mrmode=='r': results_df.to_csv('working/duplicates/dup_r_edits.csv')  
        
    return to_translate_11



def translate_targets(logger, model, to_translate, mancor):
    # attemps to rewrite reactions in a model, in a way that duplicate metabolites are solved.

    
    # define key objects: 
    modeled_mids  = [m.id for m in model.metabolites]
    results_df = []  # tabular results
    rs_to_remove = []  # to limit stoichiometric inconsistencies
    
    
    # for each reaction in the model, solve each eventual duplicate
    for r in model.reactions:
        
        # get mids, puremids, compartements, for each involved metabolite. 
        r_mids     = [m.id                   for m in r.metabolites]
        r_puremids = [m.id.rsplit('_', 1)[0] for m in r.metabolites]
        r_comps  =   [m.id.rsplit('_', 1)[1] for m in r.metabolites]
        
        
        # get number of metabolite IDs involved (compartment included). 
        # It will later used to check the presence of the exact same metabolite to both the sides of the reaction.
        # In these cases, the reaction will be removed to prevent stoichiometric inconsistency. 
        rids_involved_old = set([m.id for m in r.metabolites])
        
        
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
            
            
            # apply manual corrections if any 
            if r.id in mancor['reactions'].keys(): reaction_new = mancor['reactions'][r.id]
            
            
            # apply updated reaction
            r.build_reaction_from_string(reaction_new)
            
            
            # compute again the metabolites involved to prevent stoichiometric inconsistencies
            to_remove = '-'
            rids_involved_new = set([m.id for m in r.metabolites])
            if len(rids_involved_new) != len(rids_involved_old):  # number of metabolite IDs should not change. 
                rs_to_remove.append(r)
                to_remove = 'removed'


            # populate tabular results
            translated_mids =  f'{len(translated_mids)}: {translated_mids}'
            added_mids = f'{len(added_mids)}: {added_mids}'
            if added_mids == '0: []': added_mids = '-'
            bal_suggestions = str(r.check_mass_balance())
            if bal_suggestions == '{}': bal_suggestions = '-'
            if r.id.startswith('EX_') and len(r.metabolites)==1: bal_suggestions = '-'
            results_df.append({
                'rid': r.id, 'reaction_old': reaction_old, '~~>': '~~>', 'reaction_new': reaction_new, 
                'translated_mids': translated_mids, 'added_mids': added_mids, 'bal_suggestions': bal_suggestions, 'to_remove': to_remove})
    
    
    # remove reactions if needed
    model.remove_reactions(rs_to_remove)
    
    
    # save tabular results: 
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/duplicates/dup_m_translations.csv')



def parse_disconnected_metabolites(logger, model):
    
    
    # remove metabolites involved to 0 reactions. 
    to_remove = []
    for m in model.metabolites:
        if len(m.reactions) == 0: 
            to_remove.append(m)
    model.remove_metabolites(to_remove)
    logger.debug(f"Removed {len(to_remove)} disconnected metabolites: {[m.id for m in to_remove]}.")



def remove_duplicated_and_set_gpr(model, to_translate): 
    
    
    # first delete duplicated reactions: 
    to_delete = [model.reactions.get_by_id(rid) for rid in to_translate.keys()]
    model.remove_reactions(to_delete)
    
    
    # then update gpr where needed: 
    results_df = pnd.read_csv('working/duplicates/dup_r_edits.csv', index_col=0)
    results_df = results_df[results_df['same_geneset']==False]
    results_df = results_df.reset_index(drop=True)
    groups = results_df.groupby('replacement').groups
    for replacement in groups:
        group = results_df.iloc[groups[replacement], ]
        gprs = set(group['d_gpr'].to_list() + group['r_gpr'].to_list())
        final_gpr = '(' + ') or ('.join(list(gprs) ) + ')'
        r = model.reactions.get_by_id(replacement)
        r.gene_reaction_rule = final_gpr
        r.update_genes_from_gpr() 
                                                  


def solve_duplicates(logger, identity, coverage, refmodel, mancor_filepath):
    
    
    # log some message: 
    logger.info("Detecting duplicates using MetaNetX annotations...")
    
    
    # create subdirs without overwriting: 
    os.makedirs('working/duplicates/', exist_ok=True)
    
    
    # check presence of already computed files 
    if os.path.exists(f'working/duplicates/draft_panmodel.json') and os.path.exists(f'working/duplicates/md5.pickle'):
        if os.path.exists(f'working/duplicates/draft_panmodel_da.json') and os.path.exists(f'working/duplicates/md5_da.pickle'):
            if os.path.exists(f'working/duplicates/draft_panmodel_da_dd.json') and os.path.exists(f'working/duplicates/md5_da_dd.pickle'):
                with open('working/duplicates/md5.pickle', 'rb') as handler:
                    md5 = pickle.load(handler)
                with open('working/duplicates/md5_da.pickle', 'rb') as handler:
                    md5_da = pickle.load(handler)
                with open('working/duplicates/md5_da_dd.pickle', 'rb') as handler:
                    md5_da_dd = pickle.load(handler)
                # compare md5:
                if md5 == md5_da == md5_da_dd == get_md5_string('working/duplicates/draft_panmodel.json'):
                    # log some message: 
                    logger.info('Found all the needed files already computed. Skipping this step.')
                    # signal to skip this module:
                    return 0
    
    
    # load the needed models: 
    draft_panmodel = cobra.io.load_json_model('working/duplicates/draft_panmodel_da.json')
    reffree = cobra.io.load_json_model(f'working/free/draft_panmodel_{identity}_{coverage}.json')
    if refmodel != '-':
        refmodel_basename = os.path.basename(refmodel)
        refmodel = cobra.io.load_json_model(f'working/brh/{refmodel_basename}.refmodel_translated.json')
    else: refmodel = cobra.Model('__EMPTY__')
    
    
    # get the manual corrections dictionary:
    if mancor_filepath != '-': # existence of the file was verified above
        response = mancor_to_dict(logger, mancor_filepath)  # convert the filepath to dict
        if type(response) == int:
            if response == 1: return 1
        else: # good corrections file: 
            mancor = response  
            logger.info(  # log some message
                f"Using the provided manual corrections ({len(mancor['formulas'].keys())} formula corrections, {len(mancor['charges'].keys())} charge corrections, " + \
                f"{len(mancor['reactions'].keys())} reaction corrections, and {len(mancor['blacklist'])} reactions in blacklist)...")
    else: mancor = {'formulas': {}, 'charges': {}, 'reactions': {}, 'blacklist': []}  # emtpy, easier to handle
        
    
    # solve duplicated metabolites:
    logger.info("Detecting duplicate metabolites using MetaNetX annotations...")
    to_translate = get_translation_dictionary_mnx(draft_panmodel, mrmode='m')  # 1-to-many 
    to_translate = sort_translation_dictionary(to_translate, draft_panmodel, refmodel, reffree, mrmode='m')  # 1-to-1
    translate_targets(logger, draft_panmodel, to_translate, mancor)
    parse_disconnected_metabolites(logger, draft_panmodel)
    
    
    ### TMP
    from memote.support.consistency import find_unconserved_metabolites
    print(find_unconserved_metabolites(draft_panmodel))
    ### TMP
    
    
    
    # solve duplicated reactions:
    logger.info("Detecting duplicate reactions using MetaNetX annotations...")
    to_translate = get_translation_dictionary_mnx(draft_panmodel, mrmode='r')  # 1-to-many 
    to_translate = sort_translation_dictionary(to_translate, draft_panmodel, refmodel, reffree, mrmode='r')  # 1-to-1
    remove_duplicated_and_set_gpr(draft_panmodel, to_translate)
    
    
    # save deduplicated model
    cobra.io.save_json_model(draft_panmodel, f'working/duplicates/draft_panmodel_da_dd.json')
    
        
    # trace the parent model md5:
    parent_md5 = get_md5_string('working/duplicates/draft_panmodel.json')
    with open('working/duplicates/md5_da_dd.pickle', 'wb') as handle:
        pickle.dump(parent_md5, handle)
        
        
    return 0