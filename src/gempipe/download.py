import os
import subprocess
import glob
import pickle 
import shutil 


import pandas as pnd


def get_genomes(logger, taxids, processes):
    
    
    # create a sub-directory without overwriting
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # check the presence of already availables genomes:
    found_genomes = glob.glob('working/genomes/downloaded/*.fna')
    if len(found_genomes) > 0:
        logger.info(f"Found {len(found_genomes)} genome assemblies already stored in your ./working/ directory: skipping the download from NCBI.")
        logger.debug(f"Genomes found: " + str(found_genomes))
        return 0
        

    # execute the download
    logger.info("Downloading from NCBI all the genome assemblies linked to the provided taxids...")
    with open('working/genomes/stdout_download.txt', 'w') as stdout, open('working/genomes/stderr_download.txt', 'w') as stderr: 
        command = f"""ncbi-genome-download \
            --no-cache \
            --metadata-table working/genomes/metadata.txt \
            --retries 100 --parallel 10 \
            --output-folder working/genomes/ \
            --species-taxids {taxids} \
            --formats assembly-stats,fasta \
            --section genbank \
            bacteria"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    logger.debug("Download finished. Logs are stored in ./working/genomes/stdout_download.txt and ./working/genomes/stderr_download.txt.") 
    
    
    # moving the genomes to the right directory
    os.makedirs('working/genomes/downloaded/', exist_ok=True)
    for file in glob.glob('working/genomes/genbank/bacteria/*/*.fna.gz'):
        accession = file.split('/')[-2]
        shutil.copy(file, f'working/genomes/downloaded/{accession}.fna.gz')
    shutil.rmtree('working/genomes/genbank/') # delete the old tree
    logger.debug("Moved the downloaded genomes to ./working/genomes/downloaded/.") 
    
    
    # execute the decompression
    logger.info("Decompressing the genomes using pigz...")
    with open('working/genomes/stdout_decompression.txt', 'w') as stdout, open('working/genomes/stderr_decompression.txt', 'w') as stderr: 
        command = f"""unpigz -p {processes} working/genomes/downloaded/*.fna.gz""" 
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    logger.debug("Decompression finished. Logs are stored in ./working/genomes/stdout_decompression.txt and ./working/genomes/stderr_decompression.txt.") 
    
    
    return 0 
    
    
    
def get_metadata_table(logger):
    
    
    logger.info("Creating the metadata table for your genomes...") 
    
    
    # check if raw metadata file exists:
    if not os.path.exists('working/genomes/metadata.txt'):
        logger.error("The list of genome downloaded is missing. Please run again with -o/--overwrite option.")
        return 1
    
    
    # read the raw metadata: 
    metadata = pnd.read_csv("working/genomes/metadata.txt", sep='\t')
    
    
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
    os.makedirs('working/tables/', exist_ok=True)
    metadata.to_csv("working/tables/downloaded_genomes.csv")
    logger.info("Metadata table saved in ./working/tables/downloaded_genomes.csv.") 
    
    
    # create species-to-genome dictionary:
    species_to_genome = {}
    groups = metadata.groupby('organism_name').groups
    for species in groups.keys():
        indexes = groups[species]
        subset_metadata = metadata.iloc[indexes, ]
        species_to_genome[species] = [f'working/genomes/downloaded/{accession}.fna' for accession in subset_metadata['assembly_accession']]
    logger.debug(f"Created the species-to-genome dictionary: {str(species_to_genome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome, file)
    logger.debug(f"Saved the species-to-genome dictionary to file: ./working/genomes/species_to_genome.pickle.")
    
    
    return 0
    
    
    
def handle_manual_genomes(logger, genomes):
    
    
    # create a species-to-genome dictionary
    species_to_genome = {}
    logger.debug(f"Checking the formatting of the provided -g/-genomes attribute...") 
    
    
    # check if the user specified a folder:
    if os.path.exists(genomes):
        if os.path.isdir(genomes):
            if genomes[-1] != '/': genomes = genomes + '/'
            files = glob.glob(genomes + '*')
            species_to_genome['Spp'] = files
            ciccio
    
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
    logger.info(f"Inputted {len(species_to_genome.keys())} species with well-formatted paths to assemblies.") 
    
    
    # move the genomes to the usual directory: 
    os.makedirs('working/genomes/provided/', exist_ok=True)
    for species in species_to_genome.keys():
        copied_files = []
        for file in species_to_genome[species]:
            shutil.copy(file, 'working/genomes/provided/')
            basename = os.path.basename(file)
            copied_files.append('working/genomes/provided/' + basename)
        species_to_genome[species] = copied_files
    logger.debug(f"Input genomes copied to ./working/genomes/provided/.")
    logger.debug(f"Created the species-to-genome dictionary: {str(species_to_genome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/genomes/species_to_genome.pickle', 'wb') as file:
        pickle.dump(species_to_genome, file)
    logger.debug(f"Saved the species-to-genome dictionary to file: ./working/genomes/species_to_genome.pickle.")
    
    
    return 0
    