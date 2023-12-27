import os
import subprocess
import glob
import pickle 
import shutil 


import pandas as pnd




def get_metadata_table(logger):
    
    
    # read the raw metadata: 
    logger.info("Creating the metadata table for your genomes...") 
    metadata = pnd.read_csv("working/genomes/raw_ncbi.csv", index_col=0)
    
    
    # this table has 2 rows for each genome, one if the '_assembly_stats.txt' row.
    # here we delete the *assembly_stats rows: 
    to_drop = []
    for index, row in metadata.iterrows(): 
        if row.local_filename.endswith('_assembly_stats.txt'):
            to_drop.append(index)
    metadata = metadata.drop(to_drop)
    metadata = metadata.reset_index(drop=True)
    logger.debug("Shape of the metadata table: " + str(metadata.shape))
    
    
    # merge 'infraspecific_name' and 'isolate' to a single column 'strain_isolate': 
    metadata['infraspecific_name'] = metadata['infraspecific_name'].apply(lambda x: x.replace('strain=', '') if type(x)==str and x!='na' else '')
    metadata['isolate'] = metadata['isolate'].apply(lambda x: x if type(x)==str and x!='na' else '')
    metadata['strain_isolate'] = metadata['infraspecific_name'] + metadata['isolate']
    metadata = metadata.drop(['infraspecific_name', 'isolate'], axis=1)
    
    
    # select desired columns:
    metadata = metadata[['assembly_accession', 'bioproject', 'biosample', 'excluded_from_refseq', 'refseq_category', 'relation_to_type_material', 'species_taxid', 'organism_name', 'strain_isolate', 'version_status', 'seq_rel_date', 'submitter' ]] 
    
    
    # save the metadata table to disk:
    metadata.to_csv("working/genomes/genomes.csv")
    logger.info("Metadata table saved in ./working/genomes/genomes.csv.") 
    


def create_genomes_dictionary(logger): 
    
    
    # read the metadata table
    metadata = pnd.read_csv("working/genomes/genomes.csv", index_col=0)
    
    
    # create species-to-genome dictionary:
    species_to_genome = {}
    groups = metadata.groupby('organism_name').groups
    for species in groups.keys():
        indexes = groups[species]
        subset_metadata = metadata.iloc[indexes, ]
        species_to_genome[species] = [f'working/genomes/{accession}.fna' for accession in subset_metadata['assembly_accession']]
    logger.debug(f"Created the species-to-genome dictionary: {str(species_to_genome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome, file)
    logger.debug(f"Saved the species-to-genome dictionary to file: ./working/genomes/species_to_genome.pickle.")
    


def get_genomes(logger, taxids, cores): 
    
    
    # execute the download
    logger.info("Downloading from NCBI all the genome assemblies linked to the provided taxids...")
    with open('working/logs/stdout_download.txt', 'w') as stdout, open('working/logs/stderr_download.txt', 'w') as stderr: 
        command = f"""ncbi-genome-download \
            --no-cache \
            --metadata-table working/genomes/raw_ncbi.txt \
            --retries 100 --parallel 10 \
            --output-folder working/genomes/ \
            --species-taxids {taxids} \
            --formats assembly-stats,fasta \
            --section genbank \
            bacteria"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    logger.debug("Download finished. Logs are stored in ./working/logs/stdout_download.txt and ./working/logs/stderr_download.txt.") 
    
    
    # format the metadata
    metadata = pnd.read_csv("working/genomes/raw_ncbi.txt", sep='\t')
    metadata.to_csv("working/genomes/raw_ncbi.csv")
    os.remove("working/genomes/raw_ncbi.txt")
    
    
    # moving the genomes to the right directory
    for file in glob.glob('working/genomes/genbank/bacteria/*/*.fna.gz'):
        accession = file.split('/')[-2]
        shutil.copy(file, f'working/genomes/{accession}.fna.gz')
    shutil.rmtree('working/genomes/genbank/') # delete the old tree
    logger.debug("Moved the downloaded genomes to ./working/genomes/.") 
    
    
    # execute the decompression
    logger.info("Decompressing the genomes...")
    with open('working/logs/stdout_decompression.txt', 'w') as stdout, open('working/logs/stderr_decompression.txt', 'w') as stderr: 
        command = f"""unpigz -p {cores} working/genomes/*.fna.gz""" 
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    logger.debug("Decompression finished. Logs are stored in ./working/logs/stdout_decompression.txt and ./working/logs/stderr_decompression.txt.") 
    
    
    
def check_already_downloaded():
    
    
    # get the available files: 
    found_genomes = glob.glob('working/genomes/*.fna')
    found_metadata = os.path.exists('working/genomes/raw_ncbi.csv')
    if len(found_genomes) > 0 and found_metadata:
        
        
        # get the downloaded accessions:
        accessions = []
        for genome_file in found_genomes: 
            basename = os.path.basename(genome_file)
            accession, _ = os.path.splitext(basename)
            accessions.append(accession)
            
            
        # load the metadata table:
        metadata = pnd.read_csv("working/genomes/raw_ncbi.csv", index_col=0)
        if set(accessions) == set(metadata['assembly_accession'].to_list()):
            return True
        
        
    return False
    
    

def download_genomes(logger, taxids, cores):
    
    
    # create a sub-directory without overwriting
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # check the presence of already availables genomes:
    if check_already_downloaded(): 
        found_genomes = glob.glob('working/genomes/*.fna')
        logger.info(f"Found {len(found_genomes)} genome assemblies already stored in your ./working/ directory: skipping the download from NCBI.")
        logger.debug(f"Genomes found: " + str(found_genomes))

            
        # create metadata table and genomes dictionary: 
        get_metadata_table(logger)
        create_genomes_dictionary(logger)
            
            
        return 0    
    
          
    # download from ncbi: 
    get_genomes(logger, taxids, cores)

    
    # create the metadata table and the genomes dictionary
    get_metadata_table(logger)
    create_genomes_dictionary(logger)
    
    
    return 0 
    
    

def handle_manual_genomes(logger, genomes):
    
    
    # create a sub-directory without overwriting
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # create a species-to-genome dictionary
    species_to_genome = {}
    logger.debug(f"Checking the formatting of the provided -g/-genomes attribute...") 
    
    
    # check if the user specified a folder:
    if os.path.exists(genomes):
        if os.path.isdir(genomes):
            if genomes[-1] != '/': genomes = genomes + '/'
            files = glob.glob(genomes + '*')
            species_to_genome['Spp'] = files
    
    elif '+' in genomes and '@' in genomes: 
        for species_block in genomes.split('+'):
            species, files = species_block.split('@')
            for file in files.split(','): 
                if not os.path.exists(file):
                    logger.error("The following file provided in -g/--genomes does not exists: " + file)
                    return 1
            species_to_genome[species] = files.split(',')
            
    else: # the user has just 1 species
        for file in genomes.split(','): 
            if not os.path.exists(file):
                logger.error("The following file provided in -g/--genomes does not exists: " + file)
                return 1
        species_to_genome['Spp'] = genomes.split(',')

    
    # report a summary of the parsing: 
    logger.info(f"Inputted {len(species_to_genome.keys())} species with well-formatted paths to genomes.") 
    
    
    # move the genomes to the usual directory: 
    for species in species_to_genome.keys():
        copied_files = []
        for file in species_to_genome[species]:
            basename = os.path.basename(file)
            shutil.copyfile(file, 'working/genomes/' + basename)  # just the content, not the permissions. 
            copied_files.append('working/genomes/' + basename)
        species_to_genome[species] = copied_files
    logger.debug(f"Input genomes copied to ./working/genomes/.")
    logger.debug(f"Created the species-to-genome dictionary: {str(species_to_genome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome, file)
    logger.debug(f"Saved the species-to-genome dictionary to file: ./working/genomes/species_to_genome.pickle.")
    
    
    return 0
    