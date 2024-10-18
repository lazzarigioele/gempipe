import pickle
import os
import subprocess
import multiprocessing
import itertools
import glob
import json
import shutil 


import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results
from ..commons import get_allmeta_df



def task_bmetrics(proteome, args): 
        
        
    # retrive the arguments:
    buscodb = args['buscodb']


    # get the basename without extension:
    basename = os.path.basename(proteome)
    accession, _ = os.path.splitext(basename)


    # launch the command
    with open(f'working/logs/stdout_bmetrics_{accession}.txt', 'w') as stdout, open(f'working/logs/stderr_bmetrics_{accession}.txt', 'w') as stderr: 
        command = f"""busco -f --cpu 1 --offline \
            -i {proteome} \
            --mode proteins \
            --lineage_dataset {buscodb} \
            --download_path working/bmetrics/db/ \
            --out_path working/bmetrics/ \
            --out {accession}"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()


    # return a row for the dataframe
    return [{'accession': accession, 'completed': True}]



def compute_bmetrics(logger, cores, buscodb): 
    
    
    # logger message
    logger.info("Calculating the biological metrics to filter the genomes...")
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        
    
    # check if the metrics were already computed: 
    if os.path.exists('working/filtering/bmetrics.csv'):
        bmetrics_df = pnd.read_csv('working/filtering/bmetrics.csv', index_col=0)
        presence_list = []
        for species in species_to_proteome.keys(): 
            for proteome in species_to_proteome[species]:
                # get the basename without extension:
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)
                presence_list.append(accession in bmetrics_df['accession'].to_list())
        if all(presence_list): 
            logger.info("Found all the needed files already computed. Skipping this step.")
            return 0
    
    
    # create the worlder for biological metrics: 
    os.makedirs('working/bmetrics/', exist_ok=True)
    
    
    # check if the user specified a database
    if buscodb == 'bacteria_odb10': 
        logger.warning("We strongly suggest to set a more specific Busco database instead of 'bacteria_odb10' (use -b/--buscodb). To show the available Busco databases type 'gempipe recon -b show'.")
    
    
    # assuring the presence of the specified database
    logger.debug("Downloading the specified BUSCO database...")
    with open(f'working/logs/stdout_bmetrics_dbdownload.txt', 'w') as stdout, open(f'working/logs/stderr_bmetrics_dbdownload.txt', 'w') as stderr: 
        command = f"""busco -f \
            --download_path working/bmetrics/db/ \
            --out_path working/bmetrics/ \
            --download {buscodb}"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    logger.debug(f"Download completed for {buscodb}.")
        
      
    # create items for parallelization: 
    items = []
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]: 
            items.append(proteome)
            
            
    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['accession', 'completed']), 
            itertools.repeat('accession'), 
            itertools.repeat(logger), 
            itertools.repeat(task_bmetrics),
            itertools.repeat({'buscodb': buscodb}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)  # all_df_combined can be ignored.
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 
            
    
    # logger message
    logger.debug("Finished to compute the biological metrics.")
    
    
    # collect the short summaries:
    logger.debug("Gathering the short summaries...")
    bmetrics_df = []  
    for file in glob.glob(f"working/bmetrics/*/run_{buscodb}/short_summary.json"): 
        accession = file.replace('working/bmetrics/', '')
        accession = accession.replace(f'/run_{buscodb}/short_summary.json', '')
        jsonout = json.load(open(file, 'r'))
        
        
        # handle different versions of busco
        try: C = jsonout['results']['Complete']
        except: C = jsonout['results']['Complete percentage']
        try: sC = jsonout['results']['Single copy']
        except: sC = jsonout['results']['Single copy percentage']
        try: mC = jsonout['results']['Multi copy']
        except: mC = jsonout['results']['Multi copy percentage']
        try: F = jsonout['results']['Fragmented']
        except: F = jsonout['results']['Fragmented percentage']
        try: M = jsonout['results']['Missing']
        except: M = jsonout['results']['Missing percentage']
        
        
        bmetrics_df.append({
            'accession': accession,
            'C': C,
            'Single copy': sC,
            'Multi copy': mC,
            'F': F,
            'M': M,
            'n_markers': jsonout['results']['n_markers']
        })
    bmetrics_df = pnd.DataFrame.from_records(bmetrics_df)
    os.makedirs('working/filtering/', exist_ok=True)
    bmetrics_df.to_csv('working/filtering/bmetrics.csv')
    logger.debug("Biological metrics saved to ./working/filtering/bmetrics.csv.")
    
    
    # cleaning the workspace from useless files: 
    shutil.rmtree('working/bmetrics/')
    for file in glob.glob('./busco_*.log'): os.remove(file)
    logger.debug("Removed useless files.")
    
    
    return 0



def compute_tmetrics(logger, cores):
    
    
    # logger message
    logger.info("Calculating the technical metrics to filter the genomes...")
    
    
    # load the previously created species_to_genome: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)
        
        
    # check if the metrics were already computed: 
    if os.path.exists('working/filtering/tmetrics.csv'):
        tmetrics_df = pnd.read_csv('working/filtering/tmetrics.csv', index_col=0)
        presence_list = []
        for species in species_to_genome.keys(): 
            for genome in species_to_genome[species]:
                # get the basename without extension:
                basename = os.path.basename(genome)
                accession, _ = os.path.splitext(basename)
                presence_list.append(accession in tmetrics_df['accession'].to_list())
        if all(presence_list): 
            logger.info("Found all the needed files already computed. Skipping this step.")
            return 0

        
    # create the list of genomes to evaluate: 
    genome_files = []
    for species in species_to_genome.keys(): 
        for genome in species_to_genome[species]: 
            genome_files.append(genome)

            
    # launch the command
    with open(f'working/logs/stdout_tmetrics.txt', 'w') as stdout, open(f'working/logs/stderr_tmetrics.txt', 'w') as stderr: 
        command = f"""seqkit stats \
            --tabular \
            --basename \
            --all \
            --threads {cores} \
            --out-file working/filtering/tmetrics.csv \
            {' '.join(genome_files)}"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
        
        
    # format the table:
    tmetrics_df = pnd.read_csv('working/filtering/tmetrics.csv', sep='\t')
    tmetrics_df = tmetrics_df.rename(columns={'file': 'accession', 'num_seqs': 'ncontigs'})
    tmetrics_df['accession'] = tmetrics_df['accession'].apply(lambda x: os.path.splitext(x)[0])
    tmetrics_df.to_csv('working/filtering/tmetrics.csv')
    logger.debug("Technical metrics saved to ./working/filtering/tmetrics.csv.")
                
                
    # logger message:
    logger.debug(f"Finished to compute the technical metrics.")
    
    
    return 0
    
    
    
def figure_bmetrics(logger, outdir, bad_genomes): 
    
    
    logger.info("Producing figure for biological metrics in {outdir}/figures/bmetrics.png...")
    
    df = get_allmeta_df()

    # create new col to show filtering: 
    df['excluded'] = 0
    for accession, row in df.iterrows(): 
        if accession in bad_genomes: 
            df.loc[accession, 'excluded'] = 100
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()
    
    # draw bars:
    fig, ax = plt.subplots()
    _ = sb.barplot(df, x='strain_isolate', y='C', color='C8', ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='F', color='C9', bottom=df['C'], ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='M', color='C4', bottom=df['C']+df['F'], ax=ax)
    _ = sb.barplot(df, x='strain_isolate', y='excluded', color='white', alpha=0.55, ax=ax)
    
    # set tick labels:
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legends:
    l1 = plt.legend(handles=[Patch(color=color, label=metric) for color, metric in zip(['C8','C9','C4'], ['complete','fragmented','missing'])], title='', loc='upper left', bbox_to_anchor=(1.05, 0.5))
    l2 = plt.legend(handles=[Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())], title='', loc='lower left', bbox_to_anchor=(1.05, 0.5))
    ax.add_artist(l1)  # l2 implicitly replaces l1
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    ax.set_ylabel('% of BUSCOs')
    sb.despine()
    
    
    if len(df) <= 100:
        plt.savefig(outdir + 'figures/bmetrics.png', dpi=300, bbox_inches='tight')
    else:
        logger.info("Number of genomes is >100: producing the SVG version instead {outdir}/figures/bmetrics.svg...")
        plt.savefig(outdir + 'figures/bmetrics.svg', bbox_inches='tight')

        


def figure_tmetrics(logger, outdir, bad_genomes): 
    
    
    for tmetric in ['ncontigs', 'N50', 'sum_len']: 
        logger.info(f"Producing figure for technical metric {tmetric} in {{outdir}}/figures/{tmetric}.png...")
    
    
        df = get_allmeta_df()
        df['sum_len'] = df['sum_len'].apply(lambda x: x / 1000 / 1000)  # convert to Mb
        
        # create new col to show filtering: 
        df['excluded'] = 0
        for accession, row in df.iterrows(): 
            if accession in bad_genomes: 
                df.loc[accession, 'excluded'] = df.loc[accession, tmetric]

        # define colors:
        df = df.set_index('strain_isolate', drop=False)
        colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()

        # draw bars:
        fig, ax = plt.subplots()
        _ = sb.barplot(df, x='strain_isolate', y=tmetric, hue='strain_isolate', legend=False, palette=colors, ax=ax)
        _ = sb.barplot(df, x='strain_isolate', y='excluded', color='white', alpha=0.55, ax=ax)

        # set tick labels:
        ax.tick_params(axis='x', labelrotation=90)
        [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]

        # set legend:
        plt.legend(handles=[Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())], title='', loc='center left', bbox_to_anchor=(1.05, 0.5))

        ax.figure.set_size_inches(0.2*len(df), 4)
        ax.set_ylabel(tmetric)
        sb.despine()

        
        if len(df) <= 100:
            plt.savefig(outdir + f'figures/{tmetric}.png', dpi=300, bbox_inches='tight')
        else:
            logger.info("Number of genomes is >100: producing the SVG version instead {outdir}/figures/" + f"{tmetric}.svg...")
            plt.savefig(outdir + f'figures/{tmetric}.svg', bbox_inches='tight')

            


def filter_genomes(logger, cores, buscodb, buscoM, buscoF, ncontigs, N50, outdir, nofig):
    
    
    # compoute biological metrics: 
    response = compute_bmetrics(logger, cores, buscodb)
    if response == 1: return 1
    
    
    # compute technical metrics: 
    response = compute_tmetrics(logger, cores)
    if response == 1: return 1
    
    
    # read the metrics tables
    bmetrics_df = pnd.read_csv('working/filtering/bmetrics.csv', index_col=0)
    tmetrics_df = pnd.read_csv('working/filtering/tmetrics.csv', index_col=0)
    
    
    # get the number of Busco's scingle-copy orthologs: 
    n_sco = list(set(bmetrics_df['n_markers'].to_list()))[0]
    
    
    # check the inputted buscoM / buscoF
    busco_metrics = {'buscoM': buscoM, 'buscoF': buscoF}
    for metric, value in busco_metrics.items():
        if value.endswith('%'): 
            value = value[:-1]
            value = int(value)
        else: 
            value = int(value)
            value = value / n_sco * 100
        busco_metrics[metric] = value
        
        
    # filter genomes and proteomes based on metrics
    all_genomes = set(tmetrics_df['accession'].to_list())
    bmetrics_good = set(bmetrics_df[(bmetrics_df['M'] <= busco_metrics['buscoM']) & (bmetrics_df['F'] <= busco_metrics['buscoF'])]['accession'].to_list())
    tmetrics_good = set(tmetrics_df[(tmetrics_df['ncontigs'] <= ncontigs) & (tmetrics_df['N50'] >= N50)]['accession'].to_list())
    good_genomes = bmetrics_good.intersection(tmetrics_good)
    bad_genomes = all_genomes - good_genomes
        
    
    # load the previously created dictionaries: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
    
    
    # re-writes the dictionaries:
    species_to_genome_new = {}
    species_to_proteome_new = {}
    bad_quality = []
    for species in species_to_genome.keys(): 
        species_to_genome_new[species] = []
        species_to_proteome_new[species] = []
        for genome, proteome in zip(species_to_genome[species], species_to_proteome[species]):
            basename = os.path.basename(genome)
            accession, _ = os.path.splitext(basename)
            if accession in good_genomes:
                species_to_genome_new[species].append(genome)
                species_to_proteome_new[species].append(proteome)
            else: 
                logger.debug("Found a bad quality genome: " + genome + ". Will be ignored in subsequent analysis.")
                bad_quality.append(genome)
    logger.info(f"Found {len(bad_quality)} bad quality genomes. They will be ignored in subsequent analysis. Use --verbose to see the list.")
                
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome_new, file)
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome_new, file)
        
        
    # produce plot for genome filtering 
    if not nofig: 
        os.makedirs(outdir + 'figures/', exist_ok=True)
        figure_bmetrics(logger, outdir, bad_genomes)
        figure_tmetrics(logger, outdir, bad_genomes)
    
    
    return 0
            
        
    
    