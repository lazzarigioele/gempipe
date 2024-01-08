import os
from importlib import resources
import subprocess
import cobra
import shutil
import pickle



import pandas as pnd



from ..commons import get_blast_header
from ..commons import get_retained_accessions
    


def run_bigg_aligner(logger, cores, identity, coverage):
    
    
    # some log messages
    logger.debug("PruneU: Aligning to the BiGG genes database...")
    
    
    # run the command:
    with resources.path("gempipe.assets", "bigg_proteins.dmnd") as asset_path:  # taken from CarveMe v1.5.2
        with open(f'working/logs/stdout_biggalign.txt', 'w') as stdout, open(f'working/logs/stderr_biggalign.txt', 'w') as stderr: 
            command = f"""diamond blastp --threads {cores} \
                -d {asset_path} \
                -q working/annotation/representatives.faa \
                -o working/free/alignment.tsv \
                --ultra-sensitive --top 10 --quiet \
                --outfmt 6 {get_blast_header()}"""
            process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
            process.wait()
 

    # read the alignment
    header = get_blast_header().split(' ')
    alignment = pnd.read_csv('working/free/alignment.tsv', sep='\t', names=header)


    # filter based on alignment quality: 
    logger.debug("PruneU: Filtering based on alignment quality...")
    logger.debug(f"PruneU: Using inputted thresholds: identity {identity}%, coverage {coverage}%.") 
    query = f'pident >= {identity} & qcovhsp >= {coverage} & scovhsp >= {coverage}'
    alignment_filtered = alignment.query(query)
    logger.debug(f"PruneU: HSPs reduction: {len(alignment)} -> {len(alignment_filtered)}.")
    
    
    return alignment_filtered
    
    

def get_gprs_table(logger): 
    
    
    # some log messages
    logger.debug("PruneU: Storing all possible gene combinations for each reaction in each model... ") 
    
    
    # read the table from assets: 
    with resources.path("gempipe.assets", "bigg_gprs.csv") as asset_path:  # taken from CarveMe v1.5.2
        gprs = pnd.read_csv(asset_path)
        
        
        # create a new column Model.Gene, like 'iCN718.ABAYE_RS01830'
        gprs['BiGG_gene'] = gprs.apply(lambda row: f"{row['model']}.{row['gene'][2:]}", axis=1)
        
        
        # Some observations: 
        # A) the 'protein' column must be interpreted as 'protein complex', because it can involve several components.
        # B) Each BiGG gene can appear several times because it can be involved in several protein complexes.
        # C) Each 'protein complex' can appear several times because different reactions can require the same protein complex. 
        # D1) Each reaction can appear several times because different protein complex of the same model can enable it
        # D2) Each reaction can appear several times because different models can have the same reaction .
        
        
        return gprs
    
    
    
def get_gene_scores_table(logger, alignment, gprs):
    
    
    # some log messages
    logger.debug("PruneU: Creating the 'gene scores' table... ") 
    
    
    def agglomerate_isoenzymes(annotation_group):
        
        # the output will be a single row, so take the first row and start to edit it. 
        res_row = annotation_group.iloc[0].copy() 
        
        # get equally good clusters for this bigg.gene
        query_genes = annotation_group['query_gene'].to_list()
        
        # whe just 1 cluster: 
        if len(query_genes) == 1:
            res_row['query_gene'] = query_genes[0]
            
        else: # when many alternative clusters: 
            res_row['query_gene'] = '(' + ' or '.join(sorted(query_genes)) + ')' 
            
        # take the best performing: 
        res_row['score'] = max(annotation_group['score'].to_list())
        return res_row
    
    
    # agglomerate isoenzymes (alternative clusters for the same bigg.gene):
    annotation = alignment[['qseqid', 'sseqid', 'bitscore']]
    annotation.columns = ['query_gene', 'BiGG_gene', 'score']  # to match CarveMe names
    gene2gene = annotation \
        .sort_values(by='score', ascending=False) \
        .groupby('BiGG_gene', as_index=False) \
        .apply(agglomerate_isoenzymes)

    
    # merge with the gprs table
    gene_scores = pnd.merge(gene2gene, gprs, how='right')
    
    
    # manage spontaneus genes: 
    spontaneous_score = 0.0
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT', 'G_PP_s0001'}
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'score'] = spontaneous_score
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'query_gene'] = 'spontaneous'
    
    
    return gene_scores



def get_protein_scores_table(logger, gene_scores): 
    
    
    # some log messages
    logger.debug("PruneU: Creating the 'protein scores' table... ")
    
    
    def merge_subunits(genes):

        n_subunits = len(genes)
        genes = genes.dropna()
        
        if len(genes) != n_subunits:
            return None

        else:
            protein = ' and '.join(sorted(genes))
            if len(genes) > 1:
                return '(' + protein + ')'
            else:
                return protein            

            
    def merge_subunit_scores(scores):
        return scores.fillna(0).mean()

    
    # WARNING ! One of the consequence of this fix, is that in the 'protein_scores' could 
    # appear protein_complexes with GPR == None and score != 0.
    # In the following step, the algorithm selects the max 'score' for the whole reaction
    # (see 'return scores.max()' in merge_protein_scores()).
    # So a protein_complex-reaction-model with GPR == None could have the max 'score'. 
    # This guarantees no changes in the original final reaction scores , while providing just 
    # the right GPR.

    
    # from gene to protein scores
    protein_scores = gene_scores.groupby(['protein', 'reaction', 'model'], as_index=False) \
        .agg({'query_gene': merge_subunits, 'score': merge_subunit_scores})
    
    
    # rename to form the GPR column 
    protein_scores.rename(columns={'query_gene': 'GPR'}, inplace=True)
    
    
    return protein_scores



def get_reaction_scores_table(logger, protein_scores):
    
    
    # some log messages
    logger.debug("PruneU: Creating the 'reaction scores' table... ")
    
    
    def merge_proteins(proteins):
        
        proteins = set(proteins.dropna())
        if not proteins:
            return None
        
        gpr_str = ' or '.join(sorted(proteins))
        
        if len(proteins) > 1:
            return '(' + gpr_str + ')'
        else:
            return gpr_str

    def merge_protein_scores(scores):
        return scores.max(skipna=True)

    
    # from protein to reaction scores
    reaction_scores = protein_scores.groupby(['reaction'], as_index=False) \
        .agg({'GPR': merge_proteins, 'score': merge_protein_scores}).dropna()

    
    # normalize by the median of all the scores: 
    avg_score = reaction_scores.query('score > 0')['score'].median()
    reaction_scores['normalized_score'] = (reaction_scores['score'] / avg_score).apply(lambda x: round(x, 1))
    reaction_scores = reaction_scores.sort_values(by='normalized_score', ascending=False)
    # this table does not cotain GPR == ''. 
    
    
    return reaction_scores
    
    

def get_universe_template(logger=None, staining='neg'): 
    
    
    # some log messages
    if logger != None:  # function used also outside the pipe
        logger.debug(f"PruneU: Copying the gram {staining} universe... ")
    
    
    # get a copy of the appropriate universe: 
    if staining == 'pos':  
        with resources.path("gempipe.assets", "universe_grampos.json") as asset_path:  # taken from CarveMe v1.5.2
            universe = cobra.io.load_json_model(asset_path)  
    if staining == 'neg': 
        with resources.path("gempipe.assets", "universe_gramneg.json") as asset_path:  # taken from CarveMe v1.5.2
            universe = cobra.io.load_json_model(asset_path)
            
 
    # correct compartments:
    for m in universe.metabolites:
        if   m.compartment == 'C_c':  m.compartment = 'c'
        elif m.compartment == 'C_p':  m.compartment = 'p'
        elif m.compartment == 'C_e':  m.compartment = 'e'
    universe.compartments = {'c': 'cytosol', 'p': 'periplasm', 'e': 'extracellular'}
           
        
    # some log messages
    if logger != None:  # function used also outside the pipe
        logger.debug(f"PruneU: Done, {' '.join(['G:', str(len(universe.genes)), '|', 'R:', str(len(universe.reactions)), '|', 'M:', str(len(universe.metabolites))])}.")
    
    
    return universe
    
    

def perform_universe_pruning(logger, universe, reaction_scores): 
    
    
    # some log messages
    logger.debug(f"PruneU: Keeping reactions with genetic support...")
    
    
    # get list of reactions with genetic support: 
    # Remove the 'R_' prefix (eg R_11M3ODO -> 11M3ODO)
    tabled_rids = list(reaction_scores['reaction'].str[2:])
    
    
    # edit the universe:
    to_remove = []
    for r in universe.reactions: 
        if r.id == "Growth": continue
        if r.id in tabled_rids: 
            gpr = reaction_scores[reaction_scores['reaction'] == "R_" + r.id]['GPR'].iloc[0]
            r.gene_reaction_rule = gpr
            r.update_genes_from_gpr()
        else: 
            to_remove.append(r)
    
    
    # remove_orphans: Remove orphaned genes and metabolites from the model as well (default False).
    universe.remove_reactions(to_remove, remove_orphans=True)
    
    
    # some log messages
    logger.debug(f"PruneU: Done, {' '.join(['G:', str(len(universe.genes)), '|', 'R:', str(len(universe.reactions)), '|', 'M:', str(len(universe.metabolites))])}.")
    
    
    # substracted universe:
    universe.id = 'draft_panmodel'
    return universe
    
    

def add_exchange_reactions(logger, draft_panmodel):
    
    
    # some log messages
    logger.debug(f"PruneU: Adding the exchange reactions... ")
    
    
    # Now add the exchange reactions:
    for m in draft_panmodel.metabolites: 
        
        if m.id.rsplit('_', 1)[1] == 'e': 
            
            adding = cobra.Reaction(f"EX_{m.id}")
            adding.name = f"Exchange for {m.name}"
            
            # effectively add the reaction
            draft_panmodel.add_reactions([adding])
            added = draft_panmodel.reactions.get_by_id(f"EX_{m.id}")
            
            # define the manually corrected reaction:
            r_string = f"{m.id} -->"
            added.build_reaction_from_string(r_string) 
            
            # set bounds (just for EX_reactions we specify bounds):
            added.bounds = (0, 1000)
    
    
    # some log messages:
    logger.debug(f"PruneU: Done, {' '.join(['G:', str(len(draft_panmodel.genes)), '|', 'R:', str(len(draft_panmodel.reactions)), '|', 'M:', str(len(draft_panmodel.metabolites))])}.")
    
    
    return draft_panmodel
    
    
    
def PruneU(logger, cores, staining, identity, coverage):
    
    
    # align representatives sequences on BiGG genes, and filter HSPs based on inputted thresholds: 
    alignment = run_bigg_aligner(logger, cores, identity, coverage)
    
                 
    # get the list of all gene combinations that satisfy each reaction in each model
    gprs = get_gprs_table(logger)
    
    
    # get the 'gene_scores' table:
    gene_scores = get_gene_scores_table(logger, alignment, gprs)
    
    
    # get the 'protein_scores' table: 
    protein_scores = get_protein_scores_table(logger, gene_scores)
    
    
    # get the 'reaction_scores' table: 
    reaction_scores = get_reaction_scores_table(logger, protein_scores)
    
    
    # copy the gram pos / neg universe for subtraction: 
    universe = get_universe_template(logger, staining)
    
    
    # remove reactions not present
    draft_panmodel = perform_universe_pruning(logger, universe, reaction_scores)
    
    
    # add exchange reactions
    draft_panmodel = add_exchange_reactions(logger, draft_panmodel)
    
    
    return draft_panmodel



def eggnogg_gpr_inflator(logger, panmodel): 
    
    
    # some log messages:
    logger.debug(f"PruneU: Inflating the GPRs using functional annotation...")
    
    
    # load functional annotation:
    annot = pnd.read_csv('working/annotation/pan.emapper.annotations', sep='\t', comment='#', header=None)
    annot.columns = 'query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs'.split('\t')


    # create an editable copy of the model
    panmodel_update = panmodel.copy()

    
    # create a single column with all the important annotations: 
    annot = annot.set_index('query', drop=True, verify_integrity=True)
    discriminative_fields = ['Preferred_name', 'KEGG_ko', 'KEGG_Reaction', 'EC', 'KEGG_TC']
    annot['dense'] = ''
    for col in discriminative_fields:
        annot['dense'] = annot['dense'] + annot[col] 
        
    
    # divide in groups based on the dense annotation: 
    groups = annot.groupby(['dense']).groups

    
    #  get gids modeled AND annotated: 
    modeled_gids = set([g.id for g in panmodel.genes]) - set(['spontaneous'])  # all modeled genes
    annotated_gids = set(annot.index.to_list())  # annotated genes by eggnog-mapper
    modeled_annotated_gids = modeled_gids.intersection(annotated_gids)

    
    # track the parsed genes: 
    parsed_gids = set()
    
    
    # ignore the modeled genes that had not been annotated (group == '-----'): 
    group_not_annotated = ''.join(['-' for field in discriminative_fields])
    alts_not_annotated = annot.loc[groups[group_not_annotated], ]
    for gid in alts_not_annotated.index:
        if gid in modeled_gids: 
            parsed_gids.add(gid)
            
    
    # function to create the main log strings:
    def alts_to_string(alts, discriminative_fields):
        outstring = ''
        for index in alts.index:
            outstring = outstring + '\t' + index + ': ' + str([alts.loc[index, field] for field in discriminative_fields]) + '\n'
        return outstring
    
    
    # create a log subdirectory
    os.makedirs('working/free/gpr_inflator', exist_ok=True)
    
    
    # w_handler1: gids are divided by the number of alternatives.
    # w_handler2: gids are not divided (all together). 
    w_handler2 = open(f'working/free/gpr_inflator/gid_to_cluster_all.txt', 'w') 
    
    
    # create a gid to alternatives dictionary, for later use: 
    gid_to_alts = {}
    
    
    # separate the reactions based on the number of discovered isoforms starting from 1:
    admitted_alts = 1  
    while len(modeled_annotated_gids - parsed_gids) >= 1: 

        
        # Number of admitted isoforms with the exact same eggnog-mapper annotation: 
        w_handler1 = open(f'working/free/gpr_inflator/gid_to_cluster_{admitted_alts}.txt', 'w')
        

        # iterate the modeled genes:
        cnt = 0
        for g in panmodel.genes:
            gid = g.id
            if gid == 'spontaneous':
                continue

                
            # if this modeled gene has been annotated by eggnog-mapper: 
            if gid in annotated_gids:
                
                
                # get the group 
                group = annot.loc[gid, 'dense']
                
                
                # if this gene was not annotated (group == '-----'): 
                if group == ''.join(['-' for field in discriminative_fields]):
                    continue
                    
                
                # get the alternative genes: 
                alts = annot.loc[groups[group], ]

                
                # skip if different number of alternatives: 
                if len(alts) != admitted_alts:
                    continue

                
                # track the parsed gid
                cnt += 1
                parsed_gids.add(gid)
                
                
                # write the textual logs
                print(gid, '\n', alts_to_string(alts, discriminative_fields), end='', file=w_handler1)
                print(gid, '\n', alts_to_string(alts, discriminative_fields), end='', file=w_handler2)
                
                
                # populate the gid_to_alts dict: 
                if admitted_alts > 1: 
                    gid_to_alts[gid] = set(alts.index.to_list())
                  
        
        # close the handlers and increase the number of admitted alternatives
        w_handler1.close()
        admitted_alts += 1
    w_handler2.close()

    
    # iterate each reaction and edit the GPR accordingly:
    w_handler = open(f'working/free/gpr_inflator/gpr_updates.txt', 'w')
    for r in panmodel_update.reactions: 

        
        # get the gpr: 
        gpr = r.gene_reaction_rule
        if gpr == 'spontaneous': 
            continue

            
        # force each gid to be surrounded by spaces: 
        gpr = ' ' + gpr.replace('(', ' ( ').replace(')', ' ) ') + ' '

        
        # for each gid appearing in this GPR, sobstitute it with its group of alternatives: 
        involved_gids = [g.id for g in r.genes]
        for gid in involved_gids:
            if gid in gid_to_alts.keys(): 
                gpr = gpr.replace(f' {gid} ', f' editing_{gid} ')
        for gid in involved_gids:
            if gid in gid_to_alts.keys(): 
                or_group = ' or '.join(gid_to_alts[gid])
                gpr = gpr.replace(f' editing_{gid} ', f' ({or_group}) ')

                
        # remove spaces between parenthesis
        gpr = gpr.replace(' ( ', '(').replace(' ) ', ')')
        # remove spaces at the extremes: 
        gpr = gpr[1: -1]

        
        # if the GPR changed, update the model and track the logs: 
        if gpr != r.gene_reaction_rule:
            print(r.id, file=w_handler)
            print('\toriginal:', r.gene_reaction_rule, file=w_handler)
            r.gene_reaction_rule = gpr
            # New genes are introduced. Parethesis at the extremes are removed if not necessary. 
            r.update_genes_from_gpr()
            print('\tupdated:', r.gene_reaction_rule, file=w_handler)
    w_handler.close()

    
    # some log messages:
    logger.debug(f"PruneU: Done, {' '.join(['G:', str(len(panmodel_update.genes)), '|', 'R:', str(len(panmodel_update.reactions)), '|', 'M:', str(len(panmodel_update.metabolites))])}.")
    
    
    return panmodel_update



def network_rec(logger, cores, staining, identity, coverage):
    
    
    # create subdirs without overwriting
    os.makedirs('working/free/', exist_ok=True)
    
    
    # some log messages
    logger.info('Now starting the reference-free reconstruction (codename PruneU)...')
    
    
    # work with the open source glpk_exact
    cobra_config = cobra.Configuration()
    cobra_config.solver = "glpk_exact"
    
    
    # check if the needed files are already computed
    if os.path.exists('working/free/proc_acc.pickle'):
        with open('working/free/proc_acc.pickle', 'rb') as handler:
            proc_acc = pickle.load(handler) 
        if get_retained_accessions() == proc_acc:
            if os.path.exists('working/free/alignment.tsv'):
                if os.path.exists(f'working/free/draft_panmodel_{identity}_{coverage}.json'): 
                    if os.path.exists(f'working/free/gpr_inflator/gid_to_cluster_all.txt'):
                        if os.path.exists(f'working/free/gpr_inflator/gpr_updates.txt'):
                            # log some message: 
                            logger.info('Found all the needed files already computed. Skipping this step.')
                            # signal to skip this module:
                            return 0
    
    
    # run PruneU
    draft_panmodel = PruneU(logger, cores, staining, identity, coverage)
    
    
    # run gene recovery via functional annotation reading
    draft_panmodel = eggnogg_gpr_inflator(logger, draft_panmodel) 
    
    
    # save draft pan-model to disk
    cobra.io.save_json_model(draft_panmodel, f'working/free/draft_panmodel_{identity}_{coverage}.json')
    
    
    # some log messages
    logger.info(f"Reference-free reconstruction completed: {' '.join(['G:', str(len(draft_panmodel.genes)), '|', 'R:', str(len(draft_panmodel.reactions)), '|', 'M:', str(len(draft_panmodel.metabolites))])}.")
    
    
    # make traces to keep track of the accessions processed:
    # run_bigg_aligner(), first function to be called, works on annotation/representatives.faa,
    # so copy the 'proc_acc.pickle' inside annotation/
    shutil.copyfile('working/annotation/proc_acc.pickle', 'working/free/proc_acc.pickle')
    
    
    return 0