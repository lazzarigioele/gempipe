import multiprocessing
import itertools
import os
import glob
import shutil


import pandas as pnd
import cobra


from gempipe.interface.gaps import perform_gapfilling
from gempipe.interface.gaps import import_from_universe
from gempipe.interface.gaps import get_solver


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_media_definitions
from ..commons import apply_json_medium
from ..commons import check_panmodel_growth
from ..commons import strenghten_uptakes



def task_filler(file, args):
    
    
    # get the arguments
    panmodel = args['panmodel']
    minflux = args['minflux']
    outdir = args['outdir']
    media = args['media']
    sbml = args['sbml']
    index_col = args['index_col']
    
    
    # load the strain-specific model
    model = cobra.io.load_json_model(file)

    
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
        response = apply_json_medium(model, medium)
        
    
        # first try of the gapfilling algo: 
        inserted_rids = []  # rids needed for gapfilling
        first_sol_rids = perform_gapfilling(model, panmodel, minflux=minflux, nsol=1, verbose=False)


        # if empty solution (no reactions to add: models can already grow):
        if first_sol_rids == []: 
            inserted_rids_medium[medium_name] = []
            obj_value_medium[medium_name] = model.optimize().objective_value
            status_medium[medium_name] = model.optimize().status
            

        # if a solution was found:
        elif first_sol_rids != None:
            # add the reactions: 
            for rid in first_sol_rids:
                import_from_universe(model, panmodel, rid, bounds=None, gpr='')
                inserted_rids.append(rid)
            inserted_rids_medium[medium_name] = inserted_rids
            obj_value_medium[medium_name] = model.optimize().objective_value
            status_medium[medium_name] = model.optimize().status


        else:  # if no solution was found: 
            # starting the strenghten_uptakes trick...
            # nested 'with' statement (here + gapfilling) doesn't work, so we create a dictionary to later restore edited bounds:
            exr_ori_ss = strenghten_uptakes(model)
            exr_ori_pan = strenghten_uptakes(panmodel)
            multiplier = 1


            while (first_sol_rids==None and (minflux*multiplier)<= panmodel.slim_optimize()):
                # flux trough the objective could be too low. Starting from the same 'minflux' as before,
                # we try several gapfilling at increasing flux trough the objective. Each iteration raise 1 order of magnitude.
                first_sol_rids = perform_gapfilling(model, panmodel, minflux=minflux*multiplier, nsol=1, verbose=False)
                multiplier = multiplier * 10
            # now restore the medium changes!
            for rid in exr_ori_ss.keys(): model.reactions.get_by_id(rid).lower_bound = exr_ori_ss[rid]
            for rid in exr_ori_pan.keys(): panmodel.reactions.get_by_id(rid).lower_bound = exr_ori_pan[rid]


            # if a solution was found using the strenghten_uptakes trick:
            if first_sol_rids != None:
                # add the needed reactions:
                for rid in first_sol_rids:
                    import_from_universe(model, panmodel, rid, bounds=None, gpr='')
                    inserted_rids.append(rid)
                inserted_rids_medium[medium_name] = inserted_rids
                obj_value_medium[medium_name] = model.optimize().objective_value
                status_medium[medium_name] = model.optimize().status


                # now check if this second gap-filling strategy was enough to reach 'minflux' with the given medium: 
                res = model.optimize()
                obj_value = res.objective_value
                status = res.status
                if status=='optimal' and obj_value < minflux:  # could still be 0


                    # retry the original gap-filling (a kind of polishing): 
                    first_sol_rids = perform_gapfilling(model, panmodel, minflux=minflux, nsol=1, verbose=False)
                    if first_sol_rids != None:  # additional reactions needed to satisfy the thresholds:
                        for rid in first_sol_rids:
                            import_from_universe(model, panmodel, rid, bounds=None, gpr='')
                            inserted_rids.append(rid)
                        inserted_rids_medium[medium_name] = inserted_rids
                        obj_value_medium[medium_name] = model.optimize().objective_value
                        status_medium[medium_name] = model.optimize().status


                    else:  # still no solution despite the strenghten_uptakes trick + the polishing: 
                        gapfilling_failed = medium_name
                        break
            else:  # still no solution despite the strenghten_uptakes trick:
                gapfilling_failed = medium_name
                break
                
    
    # if gapfilling failed on some medium:
    if gapfilling_failed != False: 
        return [{index_col: accession, 'R': '-', 'inserted_rids': '-', 'solver_error': f'failing gapfilling on {gapfilling_failed}', 'obj_value_gf': '-', 'status_gf': '-'}]

    
    # remove disconnected metabolites right before saving the gapfilled model: 
    to_remove = []
    for m in model.metabolites:
        if len(m.reactions) == 0: 
            to_remove.append(m)
    model.remove_metabolites(to_remove)
    
    
    # save strain specific model to disk
    n_R = len(model.reactions)
    cobra.io.save_json_model(model, f'{outdir}/{accession}.json')
    if sbml: cobra.io.write_sbml_model(model, f'{outdir}/{accession}.xml')
    
    
    # compose the new row:
    return [{index_col: accession, 'R': n_R, 'inserted_rids': inserted_rids_medium, 'solver_error': '-', 'obj_value_gf': obj_value_medium, 'status_gf': status_medium}]



def get_gapfilling_matrix(results_df, outdir, media, matrix_file='gf_matrix_strains.csv', index_col='accession'):
    
    
    # get the gapfilled rids for this specific medium:
    for medium_name, medium in media.items():
        gf_matrix = []  # list of dictionaries future dataframe
        for accession in results_df.index: 
            inserted_rids = results_df.loc[accession, 'inserted_rids']  # should be a dict
            inserted_rids = inserted_rids[medium_name]

            
            # populate the tabular results:
            if type(inserted_rids) == str:
                if inserted_rids == '-':  # model not gapfilled.
                    gf_matrix.append({index_col: accession})
            else: 
                row_dict = {}
                for rid in inserted_rids:
                    row_dict[rid] = 1
                row_dict[index_col] = accession
                gf_matrix.append(row_dict)


        # convert to dataframe: 
        gf_matrix = pnd.DataFrame.from_records(gf_matrix)
        gf_matrix = gf_matrix.set_index(index_col, drop=True)
        gf_matrix = gf_matrix.fillna(0)  # Replace missing values with 0.
        gf_matrix = gf_matrix.astype(int)  # Force from float to int.
    
    
    # save to file:
    gf_matrix.to_csv(outdir + matrix_file.replace('.csv', f'.{medium_name}.csv'))
    
            


def strain_species_filler(logger, outdir, cores, panmodel, media_filepath, minflux, sbml, level='strain'):
    
    
    # log some messages:
    logger.info(f"Gap-filling {level}-specific models...")

   
    # create output dir
    if os.path.exists(outdir + f'{level}_models_gf/'):
        # always overwriting if already existing
        shutil.rmtree(outdir + f'{level}_models_gf/')  
    os.makedirs(outdir + f'{level}_models_gf/', exist_ok=True)
    
    
    # get the list of media on which to gap-fill: 
    media = get_media_definitions(logger, media_filepath)
    if type(media)==int: return 1   # we encountered an error.


    # check if panmodel can grow on the provided media
    if check_panmodel_growth(logger, panmodel, media, minpanflux=0.001) == False:
        return 1
    

    # create items for parallelization: 
    items = []
    for file in glob.glob(outdir + f'{level}_models/*.json'):
        items.append(file)
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    if level=='strain': index_col='accession'
    else:               index_col='species'
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat([index_col] + ['R', 'solver_error']), 
            itertools.repeat(index_col), 
            itertools.repeat(logger), 
            itertools.repeat(task_filler),  # will return a new sequences dataframe (to be concat).
            itertools.repeat({'panmodel': panmodel, 'minflux': minflux, 'outdir': outdir + f'{level}_models_gf', 'media': media, 'sbml': sbml, 'index_col': index_col}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
    
    
    # join with the previous table, and save:
    if   level=='strain':  summary_file='derive_strains.csv' 
    else:                  summary_file='derive_species.csv'
    results_df = pnd.read_csv(outdir + summary_file, index_col=0)   # set accession or species as index
    if level=='species': 
        all_df_combined.index = [i.replace('_', ' ') for i in all_df_combined.index.to_list()]
    results_df = pnd.concat([results_df, all_df_combined], axis=1)
    if   level=='species':  results_df = results_df.sort_index()   # sort by species
    else:                   results_df = results_df.sort_values(by=['species', 'strain', 'niche'], ascending=True)   # sort by species
    results_df.to_csv(outdir + summary_file)   # replace summary_file.csv
    
    
    # create the gapfilling matrix starting from 'results_df'
    if   level=='strain':  get_gapfilling_matrix(results_df, outdir, media, matrix_file='gf_matrix_strains.csv', index_col=index_col)
    else:                  get_gapfilling_matrix(results_df, outdir, media, matrix_file='gf_matrix_species.csv', index_col=index_col)
    
    
    
    return 0