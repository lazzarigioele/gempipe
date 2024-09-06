import pickle
import os


import pandas as pnd
from Bio import SeqIO, SeqRecord
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


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
    
    

def get_strain_to_recovery(logger, pam_modeled):
    
    
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
    # (eg: strain_to_recovery['GCA_01'] = {'frag': {'Cluster_2', 'Cluster_3'}}
    strain_to_recovery = get_strain_to_recovery(logger, pam_modeled)
    
    
    # convert dictionary of sets to dataframe of numbers indexed by accession: 
    # plotting-dataframe generation: 
    rec_summary = [{'accession': acc} for acc, i in strain_to_recovery.items()]
    for i, acc in enumerate(strain_to_recovery.keys()): 
        for key, values in strain_to_recovery[acc].items(): rec_summary[i][key] = len(values)
    rec_summary = pnd.DataFrame.from_records(rec_summary)
    rec_summary = rec_summary.set_index('accession', drop=True, verify_integrity=True)
    
    
    logger.info("Producing figure for gene recpvery in {outdir}/figures/genes_recovered.png...")
    
    
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
    l1 = plt.legend(handles=[Patch(color=color, label=metric) for color, metric in zip(['grey','lightgrey'], ['Recovered genes','Healthy genes'])], title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    l2 = plt.legend(handles=[Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())], title='', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('modeled gene clusters')
    sb.despine()

    plt.savefig(outdir + 'figures/genes_recovered.png', dpi=300, bbox_inches='tight')



def create_recon_plots(logger, outdir):
    
    
    logger.info("Producing figures for preliminary reconstruction metrics...")
    
    # load main assets:
    draft_panmodel = read_refmodel(outdir + 'draft_panmodel.json')
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    report = pnd.read_csv(outdir + 'report.csv', index_col=0)
    
    
    # filter pam for clusters modeled in draft_panmodel
    modeled_gids = [g.id for g in draft_panmodel.genes if g.id.startswith('Cluster_')]
    pam_modeled = pam.loc[modeled_gids, ]
    
    
    figure_genes_recovered(logger, outdir, pam_modeled)
    
    
    return 0
    
    