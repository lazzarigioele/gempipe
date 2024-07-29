import os
from importlib import resources
import subprocess
import shutil


import pandas as pnd
import cobra


from .networkrec import get_gene_scores_table
from .networkrec import get_protein_scores_table
from .networkrec import get_reaction_scores_table
from .networkrec import normalize_reaction_scores
from .networkrec import get_universe_template



from ..commons import get_blast_header




def run_tcdb_aligner(logger, cores):
    
    
    # some log messages
    logger.debug("Copying the TCDB database...")
    with resources.path("gempipe.assets", "tcdb_proteins.dmnd") as asset_path:  
        shutil.copyfile(asset_path, 'working/tcdb_transporters/tcdb_proteins.dmnd')
    
    
    # some log messages
    logger.debug("Aligning to the TCDB genes database...")
    # run the command:
    with open(f'working/logs/stdout_tcdbalign.txt', 'w') as stdout, open(f'working/logs/stderr_tcdbalign.txt', 'w') as stderr: 
        command = f"""diamond blastp --threads {cores} \
            -d working/tcdb_transporters/tcdb_proteins.dmnd \
            -q working/annotation/representatives.faa \
            -o working/tcdb_transporters/tcdb_alignment.tsv \
            --ultra-sensitive --quiet \
            --outfmt 6 {get_blast_header()}"""
        # not "--top 10" to avoid competing families.
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()



def filter_alignment(logger):
    
    
    header = get_blast_header().split(' ')
    alignment = pnd.read_csv('working/tcdb_transporters/tcdb_alignment.tsv', sep='\t', names=header)
    
    identity = 45
    positivity = 60
    coverage = 80
    query = f'pident >= {identity} & ppos >= {positivity} & qcovhsp >= {coverage} & scovhsp >= {coverage}'
    alignment_filtered = alignment.query(query)
    
    logger.debug(f"HSPs reduction: {len(alignment)} -> {len(alignment_filtered)}.")
    
    return alignment_filtered



def tcdb_matching(logger, tcdb_rs, reaction_scores_normalized):
    
    
    # get transport reactions completely supported by genes:
    rids_normalized = set(list(reaction_scores_normalized['reaction']))
    logger.debug(f"Number of retrieved transporters: {len(rids_normalized)}")


    # get transport reactions with chebi annotation:
    rids_chebi = set(list(tcdb_rs[tcdb_rs['rid'].isin(rids_normalized)]['rid']))
    logger.debug(f'With a "chebi substrate annotation: {len(rids_chebi)}/{len(rids_normalized)}')
    logger.debug(f"Missing annotation for {len(rids_normalized - rids_chebi)}: {[rid[2:].replace('_','.') for rid in (rids_normalized - rids_chebi)]}")


    # compute intersection: 
    intersection = rids_chebi.intersection(rids_normalized)
    matched = tcdb_rs[tcdb_rs['rid'].isin(intersection)]
    matched = matched.reset_index(drop=True)
    matched['gpr'] = None
    for index, row in matched.iterrows(): 
        matched.loc[index, 'gpr'] = reaction_scores_normalized[reaction_scores_normalized['reaction']==row['rid']]['GPR'].iloc[0]


    return matched



def add_putative_transporters(matched, model, universe):
    
    
    matched['added'] = None
    for index, row in matched.iterrows(): 
        
        
        chebi_2_added_rstrings = {}
        
        # get mids of involved metabolites: 
        for chebi_id, rstrings in eval(row['rstrings']).items():
            
            
            # each Chebi-ID is linked to putative transport reactions (rstrings)
            # (can be more then 1 reactions, as a Chebi-ID could be associated with more then 1 BiGG-IDs).
            # With 'added_rstrings' we keep trace of which rstrings have been added to the model for each Chebi-ID. 
            # 'chebi_2_added_rstrings' will be added to each row. 
            added_rstrings = []
            
            
            for i, rstring in enumerate(rstrings): 
                mids_involved = rstring.replace(' --> ', ' + ').replace(' <=> ', ' + ')
                mids_involved = mids_involved.split(' + ')
                
                
                # check if all 'mids_involved', are already included in the model:
                missing = set()
                for mid in mids_involved: 
                    try: m = model.metabolites.get_by_id(mid)
                    except: missing.add(mid)
                    
                    
                # if at least one involved metabolite cannot be recovered from panmodel/universe,
                # then interrupt the addition of this 'rstring' setting 'to_skip=True'
                to_skip = False
                
                    
                # check if missing mids are already modeled, but in another compartment
                for mid in missing: 
                    puremid = mid.rsplit('_', 1)[0]
                    found_m = None
                    for comp in ['c', 'p', 'e']:
                        try: found_m = model.metabolites.get_by_id(puremid + f'_{comp}')
                        except: pass
                    
                    
                    # if not found, check the universe
                    if found_m == None:
                        
                        for comp in ['c', 'p', 'e']:
                            try: found_m = universe.metabolites.get_by_id(puremid + f'_{comp}')
                            except: pass
                        
                        
                        # if missing from the universe too, interrupt the addition of this 'rstring': 
                        if found_m == None:
                            to_skip = True
                            break
                            
                            
                    # the misssing metabolite was found:
                    m = cobra.Metabolite(mid)
                    m.formula = found_m.formula
                    m.charge = found_m.charge
                    m.compartment = mid.rsplit('_', 1)[-1]
                    
                    
                    # add the missing metabolite: 
                    model.add_metabolites([m])
                          
                        
                # if at least one involved metabolite cannot be recovered from panmodel/universe,
                # then interrupt the addition of this 'rstring' setting 'to_skip=True'
                if to_skip:
                    added_rstrings.append(False)
                    break
                    
                    
                # at this point, all involved metabolites are present in the model. 
                # So we can add the reaction.
                # 'rid' will be TCDB-ID ('.' to '_') + Chebi-ID (as many transporters have specificity for more then one substrate)
                rid = row['rid'] + '__' + str(chebi_id) + '__' + str(i)
                r = cobra.Reaction(rid)
                model.add_reactions([r])
                r = model.reactions.get_by_id(rid)
                
                
                # manage rstring
                r.build_reaction_from_string(rstring)
                
                
                # manage gpr
                r.gene_reaction_rule = row['gpr']
                r.update_genes_from_gpr()
                
                
                # manage bounds
                if row['rev']=='-->' : r.bounds = (0, 1000) 
                else: r.bounds = (-1000, 1000) 
                
                
                added_rstrings.append(True)
            chebi_2_added_rstrings[chebi_id] = added_rstrings
        matched.loc[index, 'added'] = str(chebi_2_added_rstrings)                                 
                                       
        
    return model, matched



def stats_tcdb(logger, matched, plot=False): 
    
    
    # number of transporters
    n_ts = len(matched)
    
    
    # number of transporters annotated with chebi 
    n_annotated = len([eval(i) for i in matched['chebi'].to_list()  if eval(i) != {}])
    
    
    # number of transporters with at least 1 bigg reaction
    n_rstrings = 0
    for index, row in matched.iterrows(): 
        increment = False
        for chebi_id, rstrings in eval(row['rstrings']).items():
            if rstrings != set():
                increment = True
        if increment: 
            n_rstrings += 1
                
                
    # number of transporters with at least 1 bigg reaction INSERTED
    n_inserted = 0
    for index, row in matched.iterrows(): 
        increment = False
        for chebi_id, added in eval(row['added']).items():
            for response in added: 
                if response == True:
                    increment = True
        if increment: 
            n_inserted += 1
            
            
    # print statistics: 
    logger.info(f"n_ts: {n_ts} | n_annotated: {n_annotated} | n_rstrings: {n_rstrings} | n_inserted: {n_inserted}")
    
    
    # make the barplot:
    if plot: 
        categories = ['n_ts', 'n_annotated', 'n_rstrings', 'n_inserted']
        values = [n_ts, n_annotated, n_rstrings, n_inserted]  
        colors = ['C0', 'C1', 'C2', 'C3']
        fig, ax = plt.subplots()
        bars = ax.bar(categories, values, color=colors)
        [ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{bar.get_height()}', ha='center', va='bottom') for bar in bars]
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'modeled tranporters: {round(n_inserted/n_ts*100, 1)}%')


        return fig
    
    
    
def add_missing_exchanges(logger, model):
    
    added_exr_rids = []
    for m in model.metabolites: 
        if m.id.endswith('_e'):
            try: r = model.reactions.get_by_id(f'EX_{m.id}')
            except: 
    
                # create the reaction
                adding = cobra.Reaction(f"EX_{m.id}")
                adding.name = f"Exchange for {m.name}"
            
                # effectively add the reaction
                model.add_reactions([adding])
                added = model.reactions.get_by_id(f"EX_{m.id}")

                # define the manually corrected reaction:
                r_string = f"{m.id} -->"
                added.build_reaction_from_string(r_string) 

                # set bounds (just for EX_reactions we specify bounds):
                added.bounds = (0, 1000)
                
                added_exr_rids.append(f"EX_{m.id}")
                
    # some log messages:
    if logger != None: logger.debug(f"Added {len(added_exr_rids)} new EXR: {added_exr_rids}")
    
    return added_exr_rids



def tcdbing_main(logger, cores, staining):
    
    
    # create subdirs without overwriting
    os.makedirs('working/tcdb_transporters/', exist_ok=True)
    
    
    # some log messages
    logger.info('Experimental feature requested: trying to build transport reactions using TCDB...')
    
    
    # check if the needed files are already computed: 
    if False:
        logger.info('Found all the needed files already computed. Skipping this step.')
        # signal to skip this module:
        return 0
        
    
    # align representatives sequences on TCDB genes:
    run_tcdb_aligner(logger, cores)
    
    
    # filter the alignment
    alignment_filtered = filter_alignment(logger)

    
    # get the 'gprm_table':
    with resources.path("gempipe.assets", "tcdb_gprs.csv") as asset_path:  
        gprm_table = pnd.read_csv(asset_path)
    gprm_table.to_csv('working/tcdb_transporters/gprm_table.csv')
    
    
    # get the 'gene_scores' table:
    gene_scores = get_gene_scores_table(logger, alignment_filtered, gprm_table)
    gene_scores.to_csv('working/tcdb_transporters/gene_scores.csv')
    
    
    # get the 'protein_scores' table: 
    protein_scores = get_protein_scores_table(logger, gene_scores)
    protein_scores.to_csv('working/tcdb_transporters/protein_scores.csv')
    
    
    # get the 'reaction_scores' table: 
    reaction_scores = get_reaction_scores_table(logger, protein_scores)
    reaction_scores.to_csv('working/tcdb_transporters/reaction_scores.csv')
    
    
    # normalize reaction scores:
    reaction_scores_normalized = normalize_reaction_scores(reaction_scores)
    reaction_scores_normalized.to_csv('working/tcdb_transporters/reaction_scores_normalized.csv')
    
    
    # matching alingment with precomputed reaction database ('tcdb_rs'):
    with resources.path("gempipe.assets", "tcdb_rs.csv") as asset_path:  
        tcdb_rs = pnd.read_csv(asset_path)
    matched = tcdb_matching(logger, tcdb_rs, reaction_scores_normalized)
    matched.to_csv('working/tcdb_transporters/matched.csv', index=False)
    
    
    # using 'matched', expand the reference draft pan-model.
    # The resulting model will REPLACE the 'working/duplicates/draft_panmodel.json'.
    # load draft panmodel:
    panmodel = cobra.io.load_json_model('working/duplicates/draft_panmodel.json')
    uni = get_universe_template(logger, staining)
    logger.debug("Adding putative transporters...")
    panmodel, matched_added = add_putative_transporters(matched, panmodel, uni)
    matched_added.to_csv('working/tcdb_transporters/matched_added.csv')
    logger.debug("Adding missing exchanges...")
    add_missing_exchanges(logger, panmodel)
    logger.debug("Replacing working/duplicates/draft_panmodel.json...")
    cobra.io.save_json_model(panmodel, 'working/duplicates/draft_panmodel.json')
                                        
    
    # print some statistics:
    stats_tcdb(logger, matched_added)
    
    
    return 0