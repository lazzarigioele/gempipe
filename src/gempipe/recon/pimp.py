import pickle
import cobra
from importlib import resources
import os


from ..commons import get_md5_string


__PIMPCACHE__ = None


def get_annotation_info(model):
    
    
    # count how many annotation database for each object type (mets, reacs):
    # for metabolites: 
    m_annotations = set()
    for m in model.metabolites: 
        for a in m.annotation.keys():
            m_annotations.add(a)
    m_len = len(m_annotations)
    m_annotations = sorted(m_annotations)
    m_annotations= ', '.join(m_annotations)
    print(f"M ({m_len}):   " + str(m_annotations))
    
    
    # for reactions: 
    r_annotations = set()
    for r in model.reactions: 
        for a in r.annotation.keys():
            r_annotations.add(a)
    r_len = len(r_annotations)
    r_annotations = sorted(r_annotations)
    r_annotations= ', '.join(r_annotations)
    print(f"R ({r_len}):   " + str(r_annotations))
    
       
    
def make_memote_compliant(something_to_others):

                
    # make small adjustment to match the Memote/Miriam requirements:           
    for bigg_id in something_to_others.keys():
        for db in something_to_others[bigg_id].keys():
            
            
            if db =='inchikey':
                newlist = []
                for i in something_to_others[bigg_id][db]:
                    newlist.append(i.replace("InChIKey=", ''))
                something_to_others[bigg_id][db] = newlist
                
                
            if db =='chebi':
                newlist = []
                for i in something_to_others[bigg_id][db]:
                    newlist.append('CHEBI:' + i)
                something_to_others[bigg_id][db] = newlist
                
                
            if db == 'hmdb':
                # HMDB codes, to be accepted by Memote (MIRIAM) , must be of 5 
                # numbers, not 7. Comapare eg the following 'HMDB0001039', 'HMDB01039'.
                newlist = []
                for i in something_to_others[bigg_id][db]:
                    if len(i.replace('HMDB', ''))== 5:
                        newlist.append(i)
                something_to_others[bigg_id][db] = newlist
                
                
            if db == 'bigg.reaction':
                # avoid annotation repetition (eg UDCPDPS and R_UDCPDPS).
                newlist = []
                for i in something_to_others[bigg_id][db]:
                    if i.startswith('R_')==False:
                        newlist.append(i)
                something_to_others[bigg_id][db] = newlist
                
                
            if db == 'ec-code':
                # Reactions can have multiple EC codes. 
                # Not allowed are EC 1) incomplete or 2) containing letters.
                # Eg here we jeep just the second: 2.5.1.M1, 2.5.1.31, 2.5.1
                newlist = []
                for i in something_to_others[bigg_id][db]:
                    if len(i.split('.'))==4 and i.replace('.', '').isnumeric():
                        newlist.append(i)
                something_to_others[bigg_id][db] = newlist
                
    
    return something_to_others
    

    
def boost_annotations(logger, x, something_to_others, mrswitch='r', overwrite=False): 
    # Take a single reaction or metabolite and improve its annotations based on prebuilt mnx dicts.
    
    
    #  Works both for metabolites and reactions.
    if mrswitch=='r':  
        mrid = x.id
        # correct for ModelSEED reaction ids: 
        if mrid.startswith('rxn') and mrid[-3]=='_':  # for example "rxn11567_c0"
            mrid = mrid[ : -3]
    else: mrid = x.id.rsplit('_', 1)[0]  # remove compartment

    
    # extract all the annots provided by mnx for this ID:
    cnt = 0
    if mrid in something_to_others.keys(): 
        full_annots = something_to_others[mrid]

        
        # format model's annotations as set: 
        if not overwrite: 
            for annot in x.annotation.keys():
                if type(x.annotation[annot]) == str: 
                    x.annotation[annot] = set([x.annotation[annot]])
                elif type(x.annotation[annot]) == list: 
                    x.annotation[annot] = set(x.annotation[annot])
                else: 
                    logger.error(f"Found strange annotation type ({str(type(x.annotation[annot]))}) for this annotation: {x.annotation[annot]}. Please contact the developer.")
                    return 1
        else: # the oser requested an overwriting
            for annot in x.annotation.keys():
                x.annotation[annot] = set()

        
        # iterate through variuos databases provided by metanetx: 
        for db in full_annots.keys():
            # no annotations provided for this database:
            if db not in x.annotation.keys():
                x.annotation[db] = set()
            # copy annotations: 
            for annot in full_annots[db]:
                if annot not in x.annotation[db]:
                    x.annotation[db].add(annot)
                    cnt += 1  # 1 new annotation added !


        # re-format annotations as lists:
        for db in x.annotation.keys():
            x.annotation[db] = list(x.annotation[db]) 
            

    return x.annotation, cnt



def pimp_my_model(logger, model, fromchilds=False, overwrite=False):
    
    
    # autodetect ID system looking for water in cytosol
    id_sys = 'bigg'
    """
    id_sys = None
    for m in model.metabolites:
        # remove cytosol annot
        mid = m.id.rsplit('_', 1)[0]
        if mid == 'h2o': 
            id_sys = 'bigg'
            break
        elif mid == 'cpd00001':
            id_sys = 'seed'
            break
    if verbose: print(f'"{id_sys}" detected.', flush=True)
    if id_sys == None: id_sys = 'bigg'
    """
    
    
    # load the prebuilt dictionary only if this was the first request: 
    global __PIMPCACHE__
    if __PIMPCACHE__ == None:
        __PIMPCACHE__ = {} # flag for 'already loaded'
        logger.debug("Loading MetaNetX prebuilt dicts... ")
        
        
        if not fromchilds:
            with resources.path("gempipe.assets", f"{id_sys}_to_others_M.pickle") as asset_path: 
                with open(asset_path, 'rb') as handler:
                    __PIMPCACHE__['something_to_others_M'] = pickle.load(handler)
            with resources.path("gempipe.assets", f"{id_sys}_to_others_R.pickle") as asset_path: 
                with open(asset_path, 'rb') as handler:
                    __PIMPCACHE__['something_to_others_R'] = pickle.load(handler)


        else:  # required 'child' annotations: 
            with resources.path("gempipe.assets", f"{id_sys}_to_others_extended_M.pickle") as asset_path: 
                with open(asset_path, 'rb') as handler:
                    __PIMPCACHE__['something_to_others_M'] = pickle.load(handler)
            with resources.path("gempipe.assets", f"{id_sys}_to_others_extended_R.pickle") as asset_path: 
                with open(asset_path, 'rb') as handler:
                    __PIMPCACHE__['something_to_others_R'] = pickle.load(handler)
                    
    
    # get what has been loaded:
    something_to_others_M = __PIMPCACHE__['something_to_others_M']
    something_to_others_R = __PIMPCACHE__['something_to_others_R']
        
    
    # apply small corrections to meet the Memote/Miriam standards:
    something_to_others_M = make_memote_compliant(something_to_others_M)
    something_to_others_R = make_memote_compliant(something_to_others_R)
                    
    
    logger.debug("Annotating metabolites...")   
    m_cnt = 0  # counter for logger
    for m in model.metabolites:
        response = boost_annotations(logger, m, something_to_others_M, 'm', overwrite)
        if response == 1: return 1
        else: m.annotation, cnt = response
        m_cnt = m_cnt + cnt
    logger.debug(f'{m_cnt} annots added.')
    
    
    logger.debug("Annotating reactions... ")
    r_cnt = 0  # counter for logger
    for r in model.reactions:
        response = boost_annotations(logger, r, something_to_others_R, 'r', overwrite)
        if response == 1: return 1
        else: r.annotation, cnt = response
        r_cnt = r_cnt + cnt
    logger.debug(f'{r_cnt} annots added.')  
    
    
    return 0
        
        

def denovo_annotation(logger):
    
    
    # log some message: 
    logger.info("Performing de-novo model annotation...")
    
    
    # set up the cobra solver
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"
    
    
    # check presence of already computed files 
    if os.path.exists(f'working/duplicates/draft_panmodel.json') and os.path.exists(f'working/duplicates/md5.pickle'):
        if os.path.exists(f'working/duplicates/draft_panmodel_da.json') and os.path.exists(f'working/duplicates/md5_da.pickle'):
            with open('working/duplicates/md5.pickle', 'rb') as handler:
                md5 = pickle.load(handler)
            with open('working/duplicates/md5_da.pickle', 'rb') as handler:
                md5_da = pickle.load(handler)
            # compare md5:
            if md5 == md5_da == get_md5_string('working/duplicates/draft_panmodel.json'):
                # log some message: 
                logger.info('Found all the needed files already computed. Skipping this step.')
                # signal to skip this module:
                return 0
    
    
    # load the final draft panmodel
    draft_panmodel = cobra.io.load_json_model('working/duplicates/draft_panmodel.json')
    
                
    # denovo annotation
    response = pimp_my_model(logger, draft_panmodel, overwrite=True)
    if response == 1: return 1

    
    # replace file on disk:
    cobra.io.save_json_model(draft_panmodel, 'working/duplicates/draft_panmodel_da.json')
    
    
    # trace the parent model md5:
    parent_md5 = get_md5_string('working/duplicates/draft_panmodel.json')
    with open('working/duplicates/md5_da.pickle', 'wb') as handle:
        pickle.dump(parent_md5, handle)
    
    
    return 0 
    
    
    