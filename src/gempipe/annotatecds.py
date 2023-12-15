import pickle
import os
import subprocess
import multiprocessing
import itertools
import shutil


import pandas as pnd


from .commons import chunkize_items
from .commons import load_the_worker
from .commons import gather_results



            
def task_annotation(genome, args):
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)


    # launch the command
    with open(f'working/logs/stdout_annot_{accession}.txt', 'w') as stdout, open(f'working/logs/stderr_annot_{accession}.txt', 'w') as stderr: 
        command = f"""prokka --force --quiet \
            --cpus 1 \
            --outdir working/proteomes/ \
            --prefix {accession} \
            --noanno \
            --norrna \
            --notrna \
            {genome}"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
        
        
    # remove useless files:
    os.remove(f'working/proteomes/{accession}.err')
    os.remove(f'working/proteomes/{accession}.ffn')
    os.remove(f'working/proteomes/{accession}.fna')
    os.remove(f'working/proteomes/{accession}.fsa')
    os.remove(f'working/proteomes/{accession}.gbk')
    os.remove(f'working/proteomes/{accession}.log')
    os.remove(f'working/proteomes/{accession}.sqn')
    os.remove(f'working/proteomes/{accession}.tbl')
    os.remove(f'working/proteomes/{accession}.tsv')
    os.remove(f'working/proteomes/{accession}.txt')
        
    
    # return a row for the dataframe
    return {'accession': accession, 'completed': True}



def create_species_to_proteome(logger):
    
    
    # load the previously created species_to_genome: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)
        
        
    # create the species-to-proteome dictionary:
    species_to_proteome = {}
    for species in species_to_genome.keys(): 
        species_to_proteome[species] = []
        for genome in species_to_genome[species]: 
            basename = os.path.basename(genome)
            accession, _ = os.path.splitext(basename)
            species_to_proteome[species].append(f'working/proteomes/{accession}.faa')
    logger.debug(f"Created the species-to-proteome dictionary: " + str(species_to_proteome))
    
            
    # save the dictionary to disk: 
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome, file)
    logger.debug(f"Saved the species-to-proteome dictionary to file: ./working/proteome/species_to_proteome.pickle.")
    



def extract_cds(logger, cores):
    
    
    # create sub-directory without overwriting:
    logger.info("Extracting the CDSs from the genomes...")
    os.makedirs('working/proteomes/', exist_ok=True)


    # load the previously created species_to_genome: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)


    # create items for parallelization: 
    items = []
    for species in species_to_genome.keys(): 
        for genome in species_to_genome[species]: 
            items.append(genome)
            
            
    # check if the corresponding proteomes are already available: 
    already_computed = []
    for genome in items: 
        basename = os.path.basename(genome)
        accession, _ = os.path.splitext(basename)
        already_computed.append(os.path.exists(f'working/proteomes/{accession}.faa'))
    if all(already_computed):
        logger.info("Found all the proteomes already stored in your ./working/ directory: skipping this step.")
        # save the species_to_proteome dictionary to disk:
        create_species_to_proteome(logger)
        return 0
    

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
            itertools.repeat(task_annotation),
            itertools.repeat({}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # save tabular results:
    all_df_combined.to_csv('working/logs/mptab_annotatecds.csv')
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join()  
    
    
    # save the species_to_proteome dictionary to disk.
    create_species_to_proteome(logger)
    
    
    return 0



def handle_manual_proteomes(logger, proteomes):
    
    
    # create a species-to-genome dictionary
    species_to_proteome = {}
    logger.debug(f"Checking the formatting of the provided -p/-proteomes attribute...") 
    
    
    # check if the user specified a folder:
    if os.path.exists(proteomes):
        if os.path.isdir(proteomes):
            if proteomes[-1] != '/': proteomes = proteomes + '/'
            files = glob.glob(proteomes + '*')
            species_to_proteome['Spp'] = files
    
    elif '+' in proteomes and '@' in proteomes: 
        for species_block in proteomes.split('+'):
            species, files = species_block.split('@')
            for file in files.split(','): 
                if not os.path.exists(file):
                    logger.error("The following file provided in -p/--proteomes does not exists: " + file)
                    return 1
            species_to_proteome[species] = files.split(',')
            
    else: # the user has just 1 species
        for file in proteomes.split(','): 
            if not os.path.exists(file):
                logger.error("The following file provided in -p/--proteomes does not exists: " + file)
                return 1
        species_to_proteome['Spp'] = proteomes.split(',')

    
    # report a summary of the parsing: 
    logger.info(f"Inputted {len(species_to_proteome.keys())} species with well-formatted paths to proteomes.") 
    
    
    # move the genomes to the usual directory: 
    os.makedirs('working/proteomes/', exist_ok=True)
    for species in species_to_proteome.keys():
        copied_files = []
        for file in species_to_proteome[species]:
            shutil.copy(file, 'working/proteomes/')
            basename = os.path.basename(file)
            copied_files.append('working/proteomes/' + basename)
        species_to_proteome[species] = copied_files
    logger.debug(f"Input proteomes copied to ./working/proteomes/.")
    logger.debug(f"Created the species-to-proteome dictionary: {str(species_to_proteome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome, file)
    logger.debug(f"Saved the species-to-proteome dictionary to file: ./working/proteomes/species_to_proteome.pickle.")
    
    
    return 0



