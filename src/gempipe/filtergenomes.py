import pickle
import os
import subprocess
import glob
import json
import shutil 


import pandas as pnd



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
    
    
    # check if the user specified a database
    if buscodb == 'bacteria_odb10': 
        logger.warning("We strongly suggest to set a more specific Busco database instead of 'bacteria_odb10' (use -b/--buscodb). To show the available Busco databases type 'gempipe recon -b show'.")
    
    
    # create the worlder for biological metrics: 
    os.makedirs('working/bmetrics/', exist_ok=True)
        
        
    # execute busco over each genome: 
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]: 
            
            
            # get the basename without extension:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            
            
            # launch the command
            with open(f'working/logs/stdout_bmetrics_{accession}.txt', 'w') as stdout, open(f'working/logs/stderr_bmetrics_{accession}.txt', 'w') as stderr: 
                command = f"""busco -f --cpu {cores} \
                    -i {proteome} \
                    --mode proteins \
                    --lineage_dataset {buscodb} \
                    --download_path working/bmetrics/db/ \
                    --out_path working/bmetrics/ \
                    --out {accession}"""
                process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
                process.wait()
                
                
            # logger message:
            logger.debug(f"Done for {accession}.")
            
    
    # logger message
    logger.debug("Finished to compute the biological metrics.")
    
    
    # collect the short summaries:
    logger.debug("Gathering the short summaries...")
    bmetrics_df = []  
    for file in glob.glob(f"working/bmetrics/*/run_{buscodb}/short_summary.json"): 
        accession = file.replace('working/bmetrics/', '')
        accession = accession.replace(f'/run_{buscodb}/short_summary.json', '')
        jsonout = json.load(open(file, 'r'))
        bmetrics_df.append({
            'accession': accession,
            'C': jsonout['results']['Complete'],
            'Single copy': jsonout['results']['Single copy'],
            'Multi copy': jsonout['results']['Multi copy'],
            'F': jsonout['results']['Fragmented'],
            'M': jsonout['results']['Missing'],
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
    



def filter_genomes(logger, cores, buscodb, buscoM, ncontigs, N50):
    
    
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
    
    
    # check the inputted buscoM
    if buscoM.endswith('%'): 
        buscoM = buscoM[:-1]
        buscoM = int(buscoM)
    else: 
        buscoM = int(buscoM)
        buscoM = buscoM / n_sco * 100

        
    # filter genomes and proteomes based on metrics
    all_genomes = set(tmetrics_df['accession'].to_list())
    bmetrics_good = set(bmetrics_df[bmetrics_df['M'] <= buscoM]['accession'].to_list())
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
                logger.info("Found a bad quality genome: " + genome + ". Will be ignored in subsequent analysis.")
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome_new, file)
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome_new, file)
    
    
    return 0
            
        
    
    