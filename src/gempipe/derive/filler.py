import multiprocessing
import itertools
import os
import glob
import shutil


import pandas as pnd
import cobra


from gempipe.curate.gaps import perform_gapfilling
from gempipe.curate.gaps import import_from_universe
from gempipe.curate.gaps import get_solver


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_media_definitions
from ..commons import apply_json_medium
from ..commons import check_panmodel_growth
from ..commons import strenghten_uptakes



def task_strainfiller(file, args):
    
    
    # get the arguments
    panmodel = args['panmodel']
    minflux = args['minflux']
    outdir = args['outdir']
    media = args['media']
    sbml = args['sbml']
    
    
    # load the strain-specific model
    ss_model = cobra.io.load_json_model(file)

    
    # get the accession
    basename = os.path.basename(file)
    accession, _ = os.path.splitext(basename)
    
    
    # define key objects:
    gapfilling_failed = False
    inserted_rids_medium = {}  # track the inserted rids for each provided medium
    obj_value_medium = {}
    status_medium = {}
    
    
    # iterate the provided media
    for medium_name, medium in media.items():


        # apply the medium recipe:
        # response should be alwys 0 as the growth of panmodel was previously assessed. 
        response = apply_json_medium(panmodel, medium)
        response = apply_json_medium(ss_model, medium)
        
    
        # first try of the gapfilling algo: 
        inserted_rids = []  # rids needed for gapfilling
        first_sol_rids = perform_gapfilling(ss_model, panmodel, minflux=minflux, nsol=1, verbose=False)


        # if empty solution (no reactions to add: models can already grow):
        if first_sol_rids == []: 
            inserted_rids_medium[medium_name] = []
            obj_value_medium[medium_name] = ss_model.optimize().objective_value
            status_medium[medium_name] = ss_model.optimize().status
            

        # if a solution was found:
        elif first_sol_rids != None:
            # add the reactions: 
            for rid in first_sol_rids:
                import_from_universe(ss_model, panmodel, rid, bounds=None, gpr='')
                inserted_rids.append(rid)
            inserted_rids_medium[medium_name] = inserted_rids
            obj_value_medium[medium_name] = ss_model.optimize().objective_value
            status_medium[medium_name] = ss_model.optimize().status


        else:  # if no solution was found: 
            # starting the strenghten_uptakes trick...
            # nested 'with' statement (here + gapfilling) doesn't work, so we create a dictionary to later restore edited bounds:
            exr_ori_ss = strenghten_uptakes(ss_model)
            exr_ori_pan = strenghten_uptakes(panmodel)
            multiplier = 1


            while (first_sol_rids==None and (minflux*multiplier)<= panmodel.slim_optimize()):
                # flux trough the objective could be too low. Starting from the same 'minflux' as before,
                # we try several gapfilling at increasing flux trough the objective. Each iteration raise 1 order of magnitude.
                first_sol_rids = perform_gapfilling(ss_model, panmodel, minflux=minflux*multiplier, nsol=1, verbose=False)
                multiplier = multiplier * 10
            # now restore the medium changes!
            for rid in exr_ori_ss.keys(): ss_model.reactions.get_by_id(rid).lower_bound = exr_ori_ss[rid]
            for rid in exr_ori_pan.keys(): panmodel.reactions.get_by_id(rid).lower_bound = exr_ori_pan[rid]


            # if a solution was found using the strenghten_uptakes trick:
            if first_sol_rids != None:
                # add the needed reactions:
                for rid in first_sol_rids:
                    import_from_universe(ss_model, panmodel, rid, bounds=None, gpr='')
                    inserted_rids.append(rid)
                inserted_rids_medium[medium_name] = inserted_rids
                obj_value_medium[medium_name] = ss_model.optimize().objective_value
                status_medium[medium_name] = ss_model.optimize().status


                # now check if this second gap-filling strategy was enough to reach 'minflux' with the given medium: 
                res = ss_model.optimize()
                obj_value = res.objective_value
                status = res.status
                if status=='optimal' and obj_value < minflux:  # could still be 0


                    # retry the original gap-filling (a kind of polishing): 
                    first_sol_rids = perform_gapfilling(ss_model, panmodel, minflux=minflux, nsol=1, verbose=False)
                    if first_sol_rids != None:  # additional reactions needed to satisfy the thresholds:
                        for rid in first_sol_rids:
                            import_from_universe(ss_model, panmodel, rid, bounds=None, gpr='')
                            inserted_rids.append(rid)
                        inserted_rids_medium[medium_name] = inserted_rids
                        obj_value_medium[medium_name] = ss_model.optimize().objective_value
                        status_medium[medium_name] = ss_model.optimize().status


                    else:  # still no solution despite the strenghten_uptakes trick + the polishing: 
                        gapfilling_failed = medium_name
                        break
            else:  # still no solution despite the strenghten_uptakes trick:
                gapfilling_failed = medium_name
                break
                
    
    # if gapfilling failed on some medium:
    if gapfilling_failed != False: 
        return [{'accession': accession, 'R': '-', 'inserted_rids': '-', 'solver_error': f'failing gapfilling on {gapfilling_failed}', 'obj_value_gf': '-', 'status_gf': '-'}]

    
    # remove disconnected metabolites right before saving the gapfilled model: 
    to_remove = []
    for m in ss_model.metabolites:
        if len(m.reactions) == 0: 
            to_remove.append(m)
    ss_model.remove_metabolites(to_remove)
    
    
    # save strain specific model to disk
    n_R = len(ss_model.reactions)
    cobra.io.save_json_model(ss_model, f'{outdir}/{accession}.json')
    if sbml: cobra.io.write_sbml_model(ss_model, f'{outdir}/{accession}.xml')
    
    
    # compose the new row:
    return [{'accession': accession, 'R': n_R, 'inserted_rids': inserted_rids_medium, 'solver_error': '-', 'obj_value_gf': obj_value_medium, 'status_gf': status_medium}]



def get_gapfilling_matrix(results_df, outdir):
    
    
    gf_matrix = []  # list of dictionaries future dataframe
    for accession in results_df.index: 
        inserted_rids = results_df.loc[accession, 'inserted_rids']
        
        
        # populate the tabular results:
        if type(inserted_rids) == str:
            if inserted_rids == '-':  # model not gapfilled.
                gf_matrix.append({'accession': accession})
        else: 
            row_dict = {}
            for rid in inserted_rids:
                row_dict[rid] = 1
            row_dict['accession'] = accession
            gf_matrix.append(row_dict)
            
    
    # convert to dataframe: 
    gf_matrix = pnd.DataFrame.from_records(gf_matrix)
    gf_matrix = gf_matrix.set_index('accession', drop=True)
    gf_matrix = gf_matrix.fillna(0)  # Replace missing values with 0.
    gf_matrix = gf_matrix.astype(int)  # Force from float to int.
    
    
    # save to file:
    gf_matrix.to_csv(outdir + 'gf_matrix.csv')
            


def strain_filler(logger, outdir, cores, panmodel, media_filepath, minflux, sbml):
    
    
    # log some messages:
    logger.info("Gap-filling strain-specific models...")

   
    # create output dir
    if os.path.exists(outdir + 'strain_models_gf/'):
        # always overwriting if already existing
        shutil.rmtree(outdir + 'strain_models_gf/')  
    os.makedirs(outdir + 'strain_models_gf/', exist_ok=True)
    
    
    # get the list of media on which to gap-fill: 
    media = get_media_definitions(logger, media_filepath)
    if type(media)==int: return 1   # we encountered an error.


    # check if panmodel can grow on the provided media
    if check_panmodel_growth(logger, panmodel, media, minpanflux=0.001) == False:
        return 1
    

    # create items for parallelization: 
    items = []
    for file in glob.glob(outdir + 'strain_models/*.json'):
        items.append(file)
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession', 'R', 'solver_error']), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_strainfiller),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'panmodel': panmodel, 'minflux': minflux, 'outdir': outdir + 'strain_models_gf', 'media': media, 'sbml': sbml}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # join with the previous table, and save:
    results_df = pnd.read_csv(outdir + 'derive_strains.csv', index_col=0)
    results_df = pnd.concat([results_df, all_df_combined], axis=1)
    results_df = results_df.sort_values(by='accession')
    results_df.to_csv(outdir + 'derive_strains.csv')
    
    
    # create the gapfilling matrix starting from 'results_df'
    get_gapfilling_matrix(results_df, outdir)
    
    
    return 0