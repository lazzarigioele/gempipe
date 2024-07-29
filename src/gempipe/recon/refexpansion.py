import os
import shutil 
import copy
import pickle


import pandas as pnd
import cobra


from ..commons import get_retained_accessions
from ..commons import get_md5_string



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
    
    
    
def check_same_fc(r, refmodel, remove_h=False):

    
    # parse each metabolite of this reaction
    comparison = []
    for m in r.metabolites: 
        
        
        # if requested, do not consider protons:
        if remove_h:
            if m.id.rsplit('_', 1)[0] == 'h':
                continue
        
        
        # check if same molecular formula and charge: 
        ref_m = refmodel.metabolites.get_by_id(m.id)
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

    
    
def search_first_synonym(r, refmodel, remove_h=False ): 
    
    
    # check if the reference contains an equal reaction based on reactants and products mids. 
    for ref_r in refmodel.reactions:
        found_synonym = check_same_mids(r, ref_r, remove_h)
        if found_synonym: 
            return ref_r

    
    # pay attention if ATPM / NGAM is returned !
    return None
    
    
    
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
    


def add_new_reaction(r, refmodel, addedms_logger, mancor, empty_gpr=False, bounds=None):
    
    
    # don't edit the model if r.id appears in blacklist:
    if r.id in mancor['blacklist']:
        return 'blacklist'  # action for the dataframe
    
    
    # Notes for manual correction of fatty-acid-ACP formulas: 
    # Take as example hexACP: 
    # C 390 H 613 O 143 N 96 P 1 S 3    (-1) <--- proposed
    # C 6 H 11 O X     (0) <--- desired 
    # C 384 H 602 O 142 N 96 P 1 S 3     <--- to subtract
    
    
    # get metabolites already modeled
    reference_mids = set([m.id for m in refmodel.metabolites])
    
    
    # check if the reference model is purely based on stoichiometry:
    pure_stoich = True
    for m in refmodel.metabolites:
        if m.charge != None or m.formula != None: 
            pure_stoich = False
            break
    
    
    # first add missing metabolites:
    for m in r.metabolites: 
        if m.id not in reference_mids: 
            
            
            # create a new metabolite copying the attributes. 
            new_m = cobra.Metabolite(m.id)
            new_m.name = m.name
            
            # define formula and charge (here some manual corrections can be made)
            pure_mid = m.id.rsplit('_', 1)[0]
            new_m.formula = m.formula if pure_mid not in mancor['formulas'].keys() else mancor['formulas'][pure_mid]
            new_m.charge = m.charge if pure_mid not in mancor['charges'].keys() else mancor['charges'][pure_mid]
            
            # if the reference model is purely based on stoichoimetry, then we 
            # do not want to ruin its Memote metrics adding new metabolites with charge/formula:
            if pure_stoich: 
                new_m.formula = None
                new_m.charge = None
            
            # actually add the metabolite
            refmodel.add_metabolites([new_m])
            new_m = refmodel.metabolites.get_by_id(m.id)
            
            # set compartment:
            new_m.compartment = new_m.id.rsplit('_', 1)[-1]   
            
            print(f'id="{new_m.id}"\tformula="{new_m.formula}"\tcharge="{new_m.charge}"\tname="{new_m.name}"', file=addedms_logger)
            
    
    # now add missing reaction: 
    # first create a new empty reaction
    new_r = cobra.Reaction(r.id)
    new_r.name = r.name
    
    # actually add the reaction
    refmodel.add_reactions([new_r])
    new_r = refmodel.reactions.get_by_id(r.id)
    
    # define the reaction (here some manual curations can be made)
    new_r_reaction = r.reaction if r.id not in mancor['reactions'].keys() else mancor['reactions'][r.id]
    new_r.build_reaction_from_string(new_r_reaction) 
    
    # set bounds
    if bounds != None: new_r.bounds = bounds
    else: new_r.bounds = r.bounds
    
    # copy GPR: 
    if not empty_gpr:
        new_r.gene_reaction_rule = r.gene_reaction_rule 
        new_r.update_genes_from_gpr()
    
    
    return 'add' # action for the dataframe 
    
    

def mancor_to_dict(logger, mancor): 
    
    
    # convert the manual corrections file to dict: 
    resdict = {'formulas': {}, 'charges': {}, 'reactions': {}, 'blacklist': []}
    with open(mancor, 'r') as file: 
        for line in file.readlines(): 
            
            line = line.strip().rstrip()
            if line == '': continue
            if line.startswith('%'): continue  # commented line
            
            
            if line.startswith('formula.'):
                line = line[len('formula.'):]
                key, value = line.split(':', 1)
                resdict['formulas'][key] = value  
                
            elif line.startswith('charge.'):
                line = line[len('charge.'):]
                key, value = line.split(':', 1)
                value = int(value)
                resdict['charges'][key] = value
                
            elif line.startswith('reaction.'):
                line = line[len('reaction.'):]
                key, value = line.split(':', 1)
                resdict['reactions'][key] = value     
                
            elif line.startswith('blacklist.'):
                key = line[len('blacklist.'):]
                resdict['blacklist'].append(key)              
                
            else: 
                logger.error(f"Not recognized key in this line of the manual corrections file (-m/--mancor): {line}.")
                return 1

    
    return resdict

 
def expand_reference(refmodel, draft_panmodel, mancor): 
    
    
    # begin the addition of new reactions to the refmodel, to form the final draft pan-model.
    reference_rids = [r.id for r in refmodel.reactions]
    results_df = []  # summarizing all the additions
    addedms_logger = open('working/expansion/added_metabolites.txt', 'w')
    for r in draft_panmodel.reactions:
        if r.id == 'Growth': continue  # using reference biomass definition
        if [g.id for g in r.genes] == ['spontaneous']: continue  # not interested in spontaneous reactions.
        gpr = r.gene_reaction_rule
        reaction = get_detailed_reaction(r)
        
        
        # define key objects: 
        action = '-'
        rid_found = '-'
        same_mids = '-'
        same_fc = '-'
        same_gids = '-'
        ref_reaction = '-'
        final_reaction = '-'
        ref_gpr = '-'
        final_gpr = '-'
        synonym = '-'
        bal_suggestions = '-'
        
        
        # check if this rid is already modeled: 
        rid_found = r.id in reference_rids
        if rid_found:
            # get the reference reaction and gpr:
            ref_r = refmodel.reactions.get_by_id(r.id)
            
            
        else: # this rid is not modeled yet
            # search for synonyms in refmodel and take the first:
            ref_r = search_first_synonym(r, refmodel, remove_h=True)
            if ref_r != None: 
                synonym = ref_r.id
        
        
        # manage the reference reaction if any:
        if ref_r != None:
            ref_gpr = ref_r.gene_reaction_rule
            ref_reaction = get_detailed_reaction(ref_r)
            final_reaction = ref_reaction
            bal_suggestions = str(ref_r.check_mass_balance())
            
            
            # check if reactants and products are the same: 
            same_mids = check_same_mids(r, ref_r, remove_h=True)
            if same_mids: 
                # check if also formula and charge correspond: 
                same_fc = check_same_fc(r, refmodel, remove_h=True)
                
                
            # check if the gpr is the same, and define final gpr: 
            same_gids = check_same_gids(r, ref_r)
            if not same_gids: 
                # if different, glue together the two gprs, and update the reference:  
                action = 'update_gpr'
                final_gpr = get_glued_gpr(r, ref_r)
                ref_r.gene_reaction_rule = final_gpr  
                ref_r.update_genes_from_gpr()
            else: # if same gene set, just keep the definition from the reference 
                action = 'ignore'
                final_gpr = ref_gpr
                             
                
        else:  # new reaction to be included:
            
            action = add_new_reaction(r, refmodel, addedms_logger, mancor)
            if action == 'add':  # not in blacklist 
                new_r = refmodel.reactions.get_by_id(r.id)
                final_reaction = get_detailed_reaction(new_r)
                bal_suggestions = str(new_r.check_mass_balance())
                
                
                """
                from memote.support.consistency import find_unconserved_metabolites 
                cobra.io.save_json_model(refmodel, '/home/jovyan/work/validation/refmodel.json')
                refmodel2 = cobra.io.load_json_model('/home/jovyan/work/validation/refmodel.json')
                print(new_r.id, find_unconserved_metabolites(refmodel), find_unconserved_metabolites(refmodel2), len(refmodel.reactions), len(refmodel2.reactions))
                """
                
         
        
        # populate the results dataframe: 
        if bal_suggestions == '{}': bal_suggestions = '-'   # put back to standard formatting
        if r.id.startswith('EX_'): bal_suggestions = '-'   # non-sense for EX_ reactions
        new_row = {
            'rid': r.id, 'action': action, 
            'reaction': reaction, 'rid_found': rid_found, 'synonym': synonym, 'ref_reaction': ref_reaction,
            'same_mids': same_mids, 'same_fc': same_fc, 'final_reaction': final_reaction, 'bal_suggestions': bal_suggestions,
            'same_gids': same_gids, 'gpr': gpr, 'ref_gpr': ref_gpr, 'final_gpr': final_gpr, }
        results_df.append(new_row)
    
    
    # save results dataframe to disk
    addedms_logger.close()
    results_df = pnd.DataFrame.from_records(results_df)
    results_df.to_csv('working/expansion/results.csv')
    return results_df

    
    
def ref_expansion(logger, refmodel, mancor, identity, coverage): 
    
    
    # log some message
    logger.info("Expanding the reference model with new reactions taken from the reference-free reconstruction...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/expansion/', exist_ok=True)
    
    
    # check the existence of the manual corrections file:
    mancor_filepath = copy.deepcopy(mancor)
    if mancor_filepath != '-': 
        if not os.path.exists(mancor_filepath): # check the input:
            logger.error(f"Provided path to the manual corrections (-m/--mancor) does not exist: {mancor_filepath}.")
            return 1
        
        
    # do the md5: 
    if mancor_filepath != '-': # existence of the file was verified above
        md5_new = get_md5_string(mancor_filepath)
    else: 
        md5_new = '-'
            
    
    # check if the output was already computed
    if os.path.exists('working/expansion/proc_acc.pickle'):
        with open('working/expansion/proc_acc.pickle', 'rb') as handler:
            proc_acc = pickle.load(handler) 
        if get_retained_accessions() == proc_acc:
            if os.path.exists(f'working/expansion/mancor.txt'):
                md5_old = get_md5_string('working/expansion/mancor.txt')
            else: md5_old = '-'  # fresh run , or previous run without mancor.
            if md5_old == md5_new:
                if os.path.exists(f'working/expansion/draft_panmodel.json'):
                    if os.path.exists(f'working/expansion/results.csv'):
                        if os.path.exists(f'working/expansion/added_metabolites.txt'):
                            logger.info('Found all the needed files already computed. Skipping this step.')
                            # signal to skip this module:
                            return 0
    
    
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

    
    # load the reference and the reference-free reconstruction.
    refmodel_basename = os.path.basename(refmodel)
    refmodel = cobra.io.load_json_model(f'working/brh/{refmodel_basename}.refmodel_translated.json')
    draft_panmodel = cobra.io.load_json_model(f'working/free/draft_panmodel_{identity}_{coverage}.json')
    
    
    # perform the main task
    results_df = expand_reference(refmodel, draft_panmodel, mancor)
    draft_panmodel_exp = refmodel  # after the expansion
    logger.info(f"Done, {' '.join(['G:', str(len(draft_panmodel_exp.genes)), '|', 'R:', str(len(draft_panmodel_exp.reactions)), '|', 'M:', str(len(draft_panmodel_exp.metabolites))])}.")
    
    
    # count how many unbalanced reactions remained
    unbalanced = results_df[results_df['rid'].str.startswith('EX_')==False]
    unbalanced = unbalanced[unbalanced['final_reaction'].str.startswith('[[UNB]] ')]
    logger.info(f"Terminated with {len(unbalanced)} unbalanced reactions (excluding EX_ reactions).")
    
    
    # finally save the new draft panmodel
    cobra.io.save_json_model(draft_panmodel_exp, 'working/expansion/draft_panmodel.json')
    logger.debug("New draft pan-model saved to 'working/expansion/draft_panmodel.json'.")
    
    
    # make traces to keep track of the accessions processed:
    # copy the one of brh/ , since the translated refmodelis the starting point
    shutil.copyfile('working/brh/proc_acc.pickle', 'working/expansion/proc_acc.pickle')
    if mancor_filepath != '-': shutil.copyfile(mancor_filepath, 'working/expansion/mancor.txt')
    else:  # delete an eventual previous version: 
        if os.path.exists(f'working/expansion/mancor.txt'):
            os.remove(f'working/expansion/mancor.txt')
    
    
    return 0
