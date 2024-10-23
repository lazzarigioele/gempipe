import pickle
import os
import multiprocessing
import itertools


import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from ..recon.reporting import get_species_to_core_rpam
from ..recon.reporting import get_overall_core_rpam
from ..recon.reporting import get_accession_to_reactions



def figure_modeled_reactions(logger, outdir):

    # get main assets: 
    rpam = pnd.read_csv(outdir + 'rpam.csv', index_col=0)
    derive_strains = pnd.read_csv(outdir + 'derive_strains.csv')
    overall_core_rpam = get_overall_core_rpam(rpam)
    species_to_core_rpam = get_species_to_core_rpam(rpam, derive_strains)
    acc_to_reactions = get_accession_to_reactions(rpam)
    
    
    # plotting-dataframe generation: 
    reaction_summary = []   # list of dicts, future df
    for accession in acc_to_reactions.keys():
        species = derive_strains.set_index('accession', drop=True, verify_integrity=True).loc[accession, 'species']
        reaction_summary.append({
            'accession': accession,
            'accessory': len(acc_to_reactions[accession] - overall_core_rpam.union(species_to_core_rpam[species])),
            'species_core': len(species_to_core_rpam[species] - overall_core_rpam),
            'overall_core': len(overall_core_rpam),
        })
    reaction_summary = pnd.DataFrame.from_records(reaction_summary).set_index('accession', drop=True, verify_integrity=True)
    
    
    logger.info("Producing figure for modeled reactions in {outdir}/figures/reactions_modeled.png...")
    
    
    # join the dataframe to get 'strain_isolate' and 'organism_name' fields.
    genomes_df = pnd.read_csv(f'working/genomes/genomes.csv', index_col=0)
    genomes_df = genomes_df.set_index('assembly_accession', drop=True, verify_integrity=True)
    # retain only quality-filtered genomes retaining the original order: 
    genomes_df = genomes_df.loc[[i for i in genomes_df.index if i in reaction_summary.index.to_list()], ]   
    df = pnd.concat([genomes_df, reaction_summary, derive_strains.set_index('accession', drop=True, verify_integrity=True)[['R', 'G']]], axis=1)
        
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()    
    
    # draw bars:
    fig, ax = plt.subplots()
    bars_core = sb.barplot(df, x='strain_isolate', y='overall_core', color='lightgray', edgecolor='white', ax=ax, )
    for i in bars_core.patches: i.set_hatch('')  # set hatch before drawing other bars
    # much more efficient to split with a for cicle:
    for number, species in enumerate(df['organism_name'].unique()):
        df_subsetted = df[df['organism_name']==species]
        _ = sb.barplot(df_subsetted, x='strain_isolate', y='species_core', color=f'C{number}', edgecolor='white', bottom=df_subsetted['overall_core'], ax=ax) 
    _ = sb.barplot(df, x='strain_isolate', y='accessory', color='gray', edgecolor='white', bottom=df['overall_core']+df['species_core'], ax=ax)
    for i, (strain, row) in enumerate(df.iterrows()):
        offset = 0.2
        ax.text(ax.get_xticks()[i]+offset, row['R'], f"   {row['G']}", ha='center', va='bottom', fontsize=8, rotation=90)
    
    # set tick labels
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legends:
    patches = [Patch(facecolor=color, label=metric, hatch=hatch, ) for color, metric, hatch in zip(['lightgrey','grey'], ['Common core','Non-core'], ['', ''])]
    l1 = plt.legend(handles=patches, title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(df['organism_name'].unique())]
    l2 = plt.legend(handles=patches, title='Species core', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('modeled reactions')
    sb.despine()

    
    os.makedirs(outdir + 'figures/', exist_ok=True)
    if len(df) <= 100:
        plt.savefig(outdir + 'figures/reactions_modeled.png', dpi=300, bbox_inches='tight')
    else:
        logger.info("Number of genomes is >100: producing the SVG version instead {outdir}/figures/reactions_modeled.svg...")
        plt.savefig(outdir + 'figures/reactions_modeled.svg', bbox_inches='tight')
        

    return 0



def create_derive_plots(logger, outdir, nofig):
    
    
    logger.info("Producing figures for reconstruction metrics...")
    
    
    # make plots:
    if not nofig: 
        figure_modeled_reactions(logger, outdir)
    
    
    return 0