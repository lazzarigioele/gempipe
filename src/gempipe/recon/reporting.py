import pickle
import os
import multiprocessing
import itertools


import pandas as pnd
from Bio import SeqIO, SeqRecord
import cobra
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import read_refmodel



def create_panmodel_proteome(logger, outdir):
    
    logger.info("Creating the reference proteome for the draft pan-model...")
    
    
    
    # A) from the panmodel
    
    # read the final draft panmodel
    draft_panmodel = read_refmodel(outdir + 'draft_panmodel.json')
    genes_to_report = set()
    for g in draft_panmodel.genes:
        if g.id == 'spontaneous':
            continue
        genes_to_report.add(g.id)
        
    
    # collect the reference sequences
    sr_list = []
    added = set()
    for record in SeqIO.parse('working/clustering/representatives.ren.faa', "fasta"):
        cluster, cds, accession = record.description.split(' ')
        if cluster in genes_to_report:
            sr = SeqRecord.SeqRecord(record.seq, id=cluster, description=f'{cds} {accession}')
            sr_list.append(sr)
            added.add(cluster)
            genes_to_report = genes_to_report - added
            
            
    # if all the sequences were recovered, write the fasta:
    if genes_to_report == set():
        with open(outdir + 'draft_panmodel.faa', 'w') as w_handler:
            count = SeqIO.write(sr_list, w_handler, "fasta")
        logger.debug(f"{len(added)} reference sequences written to " + outdir + 'draft_panmodel.faa' + '.')
        
        
        
    # B) from the PAM
    """
    # read the final draft panmodel
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    genes_to_report = set()
    for cluster, row in pam.iterrows():
        genes_to_report.add(cluster)
        
    
    # collect the reference sequences
    sr_list = []
    added = set()
    for record in SeqIO.parse('working/clustering/representatives.ren.faa', "fasta"):
        cluster, cds, accession = record.description.split(' ')
        if cluster in genes_to_report:
            sr = SeqRecord.SeqRecord(record.seq, id=cluster, description=f'{cds} {accession}')
            sr_list.append(sr)
            added.add(cluster)
            genes_to_report = genes_to_report - added
            
            
    # if all the sequences were recovered, write the fasta:
    if genes_to_report == set():
        with open(outdir + 'draft_panproteome.faa', 'w') as w_handler:
            count = SeqIO.write(sr_list, w_handler, "fasta")
        logger.debug(f"{len(added)} reference sequences written to " + outdir + 'draft_panproteome.faa' + '.')
    """
            


def create_report(logger, outdir): 
    
    report = []  # list of dicts, future dataframe
    
    
    # get the retained genomes/proteomes (post filtering):
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        for species in species_to_proteome.keys(): 
            for proteome in species_to_proteome[species]:
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)
                
                
                # populate the table: 
                report.append({'species': species, 'accession': accession})
                
                
    # save to file
    report = pnd.DataFrame.from_records(report)
    report.to_csv(outdir + 'report.csv')
    
    
    
    return 0 
    
    

def get_accession_to_recovery(logger, pam_modeled):
    
    # parse the pam (modeled clusters) to obtain metrics on gene recovery, per accession. 
    
    acc_to_recovery = {}
    for acc in pam_modeled.columns:
        acc_to_recovery[acc] = {
            'healthy': set(), 'recovered': set(),
            'frag': set(), 'refound': set(), 'overlap': set(), 
            'stop': set()}
        
        for index, row in pam_modeled.iterrows():
            cell = pam_modeled.loc[index, acc]
            if type(cell) == float: 
                continue

            one_healthy = False
            only_stop = True
            is_broken, is_refound, is_overlap = False, False, False

            for cds in cell.split(';'):
                if all(i not in cds for i in ['_stop', '_frag', '_refound', '_overlap']):
                    one_healthy = True
                if '_stop' not in cds: 
                    only_stop = False

                if '_frag' in cds: 
                    is_broken = True
                if '_refound' in cds:
                    is_refound = True
                if '_overlap' in cds: 
                    is_overlap = True

            if only_stop:    acc_to_recovery[acc]['stop'].add(index)
            else: 
                if one_healthy: acc_to_recovery[acc]['healthy'].add(index)
                else:
                    if is_broken:  acc_to_recovery[acc]['frag'].add(index)
                    if is_refound:  acc_to_recovery[acc]['refound'].add(index)
                    if is_overlap:  acc_to_recovery[acc]['overlap'].add(index)

                    if is_broken or is_refound or is_overlap:
                        acc_to_recovery[acc]['recovered'].add(index)
                        
        # could be parallelized!
        logger.debug(f"PAM parsed for {acc}.")

        
    return acc_to_recovery



def figure_genes_recovered(logger, outdir, pam_modeled):
    
    logger.info("Gathering metrics for gene recovery...")
    
    # get the dictionary of sets summarizing gene recovery 
    # (eg: acc_to_recovery['GCA_01'] = {'frag': {'Cluster_2', 'Cluster_3'}}
    acc_to_recovery = get_accession_to_recovery(logger, pam_modeled)
    
    
    # convert dictionary of sets to dataframe of numbers indexed by accession: 
    # plotting-dataframe generation: 
    rec_summary = [{'accession': acc} for acc, i in acc_to_recovery.items()]
    for i, acc in enumerate(acc_to_recovery.keys()): 
        for key, values in acc_to_recovery[acc].items(): rec_summary[i][key] = len(values)
    rec_summary = pnd.DataFrame.from_records(rec_summary)
    rec_summary = rec_summary.set_index('accession', drop=True, verify_integrity=True)
    
    
    logger.info("Producing figure for gene recovery in {outdir}/figures/genes_recovered.png...")
    
    
    # join the dataframe to get 'strain_isolate' and 'organism_name' fields.
    genomes_df = pnd.read_csv('working/genomes/genomes.csv', index_col=0)
    genomes_df = genomes_df.set_index('assembly_accession', drop=True, verify_integrity=True)
    # retain only quality-filtered genomes retaining the original order: 
    genomes_df = genomes_df.loc[[i for i in genomes_df.index if i in rec_summary.index.to_list()], ]   
    df = pnd.concat([genomes_df, rec_summary], axis=1)                          
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()
        
    # draw bars:
    fig, ax = plt.subplots()
    _ = sb.barplot(df, x='strain_isolate', y='healthy', color='lightgrey', ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='recovered', color='grey', bottom=df['healthy'], ax=ax)
    
    # set tick labels
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legends:
    l1 = plt.legend(handles=[Patch(color=color, label=metric) for color, metric in zip(['grey','lightgrey'], ['Recovered','Healthy'])], title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    l2 = plt.legend(handles=[Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())], title='', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('modeled gene clusters')
    sb.despine()

    try: 
        plt.savefig(outdir + 'figures/genes_recovered.png', dpi=300, bbox_inches='tight')
    except:  # the png image could be too large, so we produce svg
        logger.info(f"PNG was too large: producing the SVG version instead {outdir}/figures/genes_recovered.svg...")
        plt.savefig(outdir + 'figures/genes_recovered.svg', bbox_inches='tight')


    
def get_species_to_core(pam_modeled, report):
    
    # get core genes of each species (dict of sets)
    species_to_core = {}
    
    for species in report['species'].unique(): 
        species_to_core[species] = set()  # row indexes
        
        accs = report[report['species']==species]['accession'].to_list()

        for index, row in pam_modeled[accs].iterrows():
            is_core = True
            for acc in accs:   # iterate over accessions for this species
                cell = row[acc]
                if type(cell) == float: 
                    is_core = False
                    break
                else:
                    only_stop = True
                    for cds in cell.split(';'):
                        if '_stop' not in cds:
                            only_stop = False
                    if only_stop:
                        is_core = False
                        break
            if is_core:
                species_to_core[species].add(index)
                
    return species_to_core
    
    

def get_overall_core(pam_modeled):
    
    # get core genes (considering all quality-filtered genomes in input)
    overall_core = set()
    
    for index, row in pam_modeled.iterrows():
        is_core = True
        for acc in pam_modeled.columns:   # iterate over accessions for this species
            cell = row[acc]
            if type(cell)==float: 
                is_core = False
            else:
                only_stop = True
                for cds in cell.split(';'):
                    if '_stop' not in cds:
                        only_stop = False
                if only_stop:
                    is_core = False
        if is_core:
            overall_core.add(index)
            
    return overall_core



def get_accession_to_genes(pam_modeled):

    # get all gene clusters of each quality-filtered genome (dict of sets): 
    strain_to_genes = {}
    
    for acc in pam_modeled.columns:
        strain_to_genes[acc] = set()
        for index, row in pam_modeled.iterrows(): 
            cell = pam_modeled.loc[index, acc]  
            
            if type(cell) != float:
                only_stop = True
                for cds in cell.split(';'):
                    if '_stop' not in cds:
                        only_stop = False
                if not only_stop:
                    strain_to_genes[acc].add(index)
                    
    return strain_to_genes
    
    
    
def figure_modeled_genes(logger, outdir, pam_modeled, report, draft_panmodel):
    
    logger.info("Gathering metrics for modeled gene clusters...")
    
    
    # get main assets:
    overall_core = get_overall_core(pam_modeled)
    species_to_core = get_species_to_core(pam_modeled, report)
    acc_to_genes = get_accession_to_genes(pam_modeled)
    

    # plotting-dataframe generation: 
    gene_summary = []   # list of dicts, future df
    for accession in acc_to_genes.keys():
        species = report.set_index('accession', drop=True, verify_integrity=True).loc[accession, 'species']
        gene_summary.append({
            'accession': accession,
            'accessory': len(acc_to_genes[accession] - overall_core.union(species_to_core[species])),
            'species_core': len(species_to_core[species] - overall_core),
            'overall_core': len(overall_core),
        })
    gene_summary = pnd.DataFrame.from_records(gene_summary).set_index('accession', drop=True, verify_integrity=True)
    
    
    logger.info("Producing figure for modeled gene clusters in {outdir}/figures/genes_modeled.png...")
    
    
    # join the dataframe to get 'strain_isolate' and 'organism_name' fields.
    genomes_df = pnd.read_csv('working/genomes/genomes.csv', index_col=0)
    genomes_df = genomes_df.set_index('assembly_accession', drop=True, verify_integrity=True)
    # retain only quality-filtered genomes retaining the original order: 
    genomes_df = genomes_df.loc[[i for i in genomes_df.index if i in gene_summary.index.to_list()], ]   
    df = pnd.concat([genomes_df, gene_summary], axis=1)
        
    
    # handle the draft pan-GSMM in a dedicated dataframe: 
    df_pan = df.copy()
    df_pan.loc['draft pan-GSMM'] = 0  # row
    df_pan['pan'] = 0   # column
    df_pan.loc['draft pan-GSMM', 'pan'] = len(draft_panmodel.genes)   # add information about the draft pan-GSMM in a dedicated column
    df_pan.loc['draft pan-GSMM', 'strain_isolate'] = 'draft pan-GSMM'
        
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()
    colors['draft pan-GSMM'] = 'black'
    
    # draw bars:
    fig, ax = plt.subplots()
    _ = sb.barplot(df, x='strain_isolate', y='overall_core', color='grey', ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='species_core', palette=colors, bottom=df['overall_core'], hue='strain_isolate', legend=False, ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='accessory', color='lightgrey', bottom=df['overall_core']+df['species_core'], ax=ax)
    _ = sb.barplot(df_pan, x='strain_isolate', y='pan', color='black', ax=ax)
    
    # set tick labels
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legends:
    l1 = plt.legend(handles=[Patch(color=color, label=metric) for color, metric in zip(['grey','lightgrey'], ['Common core','Non-core'])], title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    handles = [Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())]
    handles.append(Patch(color=f'black', label='draft pan-GSMM'))
    l2 = plt.legend(handles=handles, title='', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('modeled gene clusters')
    sb.despine()

    try: 
        plt.savefig(outdir + 'figures/genes_modeled.png', dpi=300, bbox_inches='tight')
    except:  # the png image could be too large, so we produce svg
        logger.info(f"PNG was too large: producing the SVG version instead {outdir}/figures/genes_modeled.svg...")
        plt.savefig(outdir + 'figures/genes_modeled.svg', bbox_inches='tight')
        
    

def task_rpam_parse_column(accession, args): 
        
    # retrive the arguments:
    draft_panmodel = args['draft_panmodel']
    pam_modeled = args['pam_modeled']


    # create a copy of the model
    ss_model = draft_panmodel.copy()

    # define which genes to remove: 
    to_remove = []
    for cluster, row in pam_modeled.iterrows():
        cell = pam_modeled.loc[cluster, accession]
        if type(cell)==float: 
            to_remove.append(ss_model.genes.get_by_id(cluster))
            continue
        only_stop = True
        for cds in cell.split(';'):
            if '_stop' not in cds:
                only_stop = False
        if only_stop:
            to_remove.append(ss_model.genes.get_by_id(cluster))
            continue
    cobra.manipulation.delete.remove_genes(ss_model, to_remove, remove_reactions=True)


    # create new row:
    row_dict = {'accession': accession}
    for r in ss_model.reactions: 
        row_dict[r.id] = 1


    # return the new dataframe row
    return [row_dict]
    
    
    
def get_rpam_from_panmodel(logger, cores, draft_panmodel, pam_modeled):
    # create a reaction presence/absence matrix starting from a PAM and a draft panmodel: 
    
    
    # create items for parallelization: 
    items = []
    for accession in pam_modeled.columns: 
        items.append(accession)
        
        
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession'] + [r.id for r in draft_panmodel.reactions]), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_rpam_parse_column),
            itertools.repeat({'draft_panmodel': draft_panmodel, 'pam_modeled': pam_modeled}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)  # all_df_combined can be ignored.
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join()  
    
    
    # format the rpam
    rpam = all_df_combined
    # now replace missing values with 0. 
    # In 'older' pandas versions, fillna() was doing silent dowcasting (see https://medium.com/@felipecaballero/deciphering-the-cryptic-futurewarning-for-fillna-in-pandas-2-01deb4e411a1)
    # In 'future' pandas versions, silent downcasting will be forbidden. 
    # In 'middle' pandas versions, there is an option to simulate the future behaviour.
    # What follows should guarantee the old behaviour in 'older' and 'middle' pandas versions:
    try:
        with pnd.option_context('future.no_silent_downcasting', True):
            rpam = rpam.fillna(0).infer_objects()   # infer_objects() will do downcasting like in 'older' pandas versions
    except: # OptionError: "No such keys(s): 'future.no_silent_downcasting'" ==> 'older' pandas version
        rpam = rpam.fillna(0)  
    rpam = rpam.astype(int)  # force from float to int.
    rpam = rpam.T  # transpose: reactions as rows.
        
        
    return rpam
    


def get_species_to_core_rpam(rpam, report):
    # get core genes of each species (dict of sets)
    
    species_to_core = {}
    for species in report['species'].unique(): 
        species_to_core[species] = set()  # row indexes
        
        accs = report[report['species']==species]['accession'].to_list()

        for index, row in rpam[accs].iterrows():
            is_core = True
            for acc in accs:   # iterate over accessions for this species
                cell = row[acc]
                if cell == 0: 
                    is_core = False
                    break
            if is_core:
                species_to_core[species].add(index)
                
    return species_to_core


    
def get_overall_core_rpam(rpam):
    # given a rpam, get the core reactions (set): 
    
    
    overall_core = set()
    for index, row in rpam.iterrows():
        is_core = True
        for acc in rpam.columns:  
            cell = row[acc]
            if cell == 0: 
                is_core = False
        if is_core: 
            overall_core.add(index)
            
    return overall_core
    
    

def get_accession_to_reactions(rpam):
    # get reactions fore each quality-filtered input genome: 
    
    acc_to_reactions = {}
    for acc in rpam.columns:
        acc_to_reactions[acc] = set()
        for index, row in rpam.iterrows(): 
            cell = rpam.loc[index, acc]   # get columns (accs) for this species.
            if cell != 0:
                acc_to_reactions[acc].add(index)
                    
    return acc_to_reactions
    
    
    
def figure_modeled_reactions(logger, outdir, cores, pam_modeled, report, draft_panmodel):
    
    
    # get main assets: 
    rpam = get_rpam_from_panmodel(logger, cores, draft_panmodel, pam_modeled)
    overall_core_rpam = get_overall_core_rpam(rpam)
    species_to_core_rpam = get_species_to_core_rpam(rpam, report)
    acc_to_reactions = get_accession_to_reactions(rpam)

    
    # plotting-dataframe generation: 
    reaction_summary = []   # list of dicts, future df
    for accession in acc_to_reactions.keys():
        species = report.set_index('accession', drop=True, verify_integrity=True).loc[accession, 'species']
        reaction_summary.append({
            'accession': accession,
            'accessory': len(acc_to_reactions[accession] - overall_core_rpam.union(species_to_core_rpam[species])),
            'species_core': len(species_to_core_rpam[species] - overall_core_rpam),
            'overall_core': len(overall_core_rpam),
        })
    reaction_summary = pnd.DataFrame.from_records(reaction_summary).set_index('accession', drop=True, verify_integrity=True)
    
    
    logger.info("Producing figure for modeled reactions in {outdir}/figures/reactions_modeled.png...")
    
    
    # join the dataframe to get 'strain_isolate' and 'organism_name' fields.
    genomes_df = pnd.read_csv('working/genomes/genomes.csv', index_col=0)
    genomes_df = genomes_df.set_index('assembly_accession', drop=True, verify_integrity=True)
    # retain only quality-filtered genomes retaining the original order: 
    genomes_df = genomes_df.loc[[i for i in genomes_df.index if i in reaction_summary.index.to_list()], ]   
    df = pnd.concat([genomes_df, reaction_summary], axis=1)
        
    
    # handle the draft pan-GSMM in a dedicated dataframe: 
    df_pan = df.copy()
    df_pan.loc['draft pan-GSMM'] = 0  # row
    df_pan['pan'] = 0   # column
    df_pan.loc['draft pan-GSMM', 'pan'] = len(draft_panmodel.reactions)   # add information about the draft pan-GSMM in a dedicated column
    df_pan.loc['draft pan-GSMM', 'strain_isolate'] = 'draft pan-GSMM'
        
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()
    colors['draft pan-GSMM'] = 'black'
    
    # draw bars:
    fig, ax = plt.subplots()
    _ = sb.barplot(df, x='strain_isolate', y='overall_core', color='grey', ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='species_core', palette=colors, bottom=df['overall_core'], hue='strain_isolate', legend=False, ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='accessory', color='lightgrey', bottom=df['overall_core']+df['species_core'], ax=ax)
    _ = sb.barplot(df_pan, x='strain_isolate', y='pan', color='black', ax=ax)
    
    # set tick labels
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legends:
    l1 = plt.legend(handles=[Patch(color=color, label=metric) for color, metric in zip(['grey','lightgrey'], ['Common core','Non-core'])], title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    handles = [Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())]
    handles.append(Patch(color=f'black', label='draft pan-GSMM'))
    l2 = plt.legend(handles=handles, title='', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('modeled reactions')
    sb.despine()

    try: 
        plt.savefig(outdir + 'figures/reactions_modeled.png', dpi=300, bbox_inches='tight')
    except:  # the png image could be too large, so we produce svg
        logger.info(f"PNG was too large: producing the SVG version instead {outdir}/figures/reactions_modeled.svg...")
        plt.savefig(outdir + 'figures/reactions_modeled.svg', bbox_inches='tight')
            

    
def create_recon_plots(logger, outdir, cores, nofig):
    
    
    logger.info("Producing figures for preliminary reconstruction metrics...")
    
    # load main assets:
    draft_panmodel = read_refmodel(outdir + 'draft_panmodel.json')
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    report = pnd.read_csv(outdir + 'report.csv', index_col=0)
    
    
    # filter pam for clusters modeled in draft_panmodel
    modeled_gids = [g.id for g in draft_panmodel.genes if g.id.startswith('Cluster_')]
    pam_modeled = pam.loc[modeled_gids, ]
    
    
    # make 3 plots:
    if not nofig: 
        figure_genes_recovered(logger, outdir, pam_modeled)
        figure_modeled_genes(logger, outdir, pam_modeled, report, draft_panmodel)
        figure_modeled_reactions(logger, outdir, cores, pam_modeled, report, draft_panmodel)
    
    
    return 0
    
    