import os
import subprocess
import glob
import pickle 
import shutil 


import pandas as pnd


from ..commons import get_genomes_csv
from ..commons import remove_duplicated_strain_ids
from ..commons import update_metadata_manual



def get_metadata_table(logger, rawmeta_filepath):
    
    
    # read the raw metadata: 
    logger.info("Creating the metadata table for your genomes...") 
    metadata = pnd.read_csv(rawmeta_filepath, index_col=0)
    
    
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
    # Warning: the same columns are used in handle_manual_genomes()
    metadata = metadata[['assembly_accession', 'bioproject', 'biosample', 'excluded_from_refseq', 'refseq_category', 'relation_to_type_material', 'species_taxid', 'organism_name', 'strain_isolate', 'version_status', 'seq_rel_date', 'submitter' ]] 
    
    
    # save the metadata table to disk:
    metadata = remove_duplicated_strain_ids(metadata)
    metadata = metadata.sort_values(by=['organism_name', 'strain_isolate'], ascending=True)   # sort by species
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
    
    
    # get metadata basename: 
    taxids_sorted = sorted(taxids.split(','))
    meta_basename = f"raw_ncbi_{'_'.join(taxids_sorted)}"
    
    
    # execute the download
    logger.info("Downloading from NCBI all the genome assemblies linked to the provided taxids...")
    with open('working/logs/stdout_download.txt', 'w') as stdout, open('working/logs/stderr_download.txt', 'w') as stderr: 
        command = f"""ncbi-genome-download \
            --no-cache \
            --metadata-table working/genomes/{meta_basename}.txt \
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
    metadata = pnd.read_csv(f"working/genomes/{meta_basename}.txt", sep='\t')
    metadata.to_csv(f"working/genomes/{meta_basename}.csv")
    os.remove(f"working/genomes/{meta_basename}.txt")
    
    
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
    
    

def download_genomes(logger, taxids, cores, metadata_man):
    
    
    # create a sub-directory without overwriting
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # get metadata basename: 
    taxids_sorted = sorted(taxids.split(','))
    meta_basename = f"raw_ncbi_{'_'.join(taxids_sorted)}"
    
    
    # check if the genomes were already downloaded:
    if os.path.exists(f'working/genomes/{meta_basename}.csv'):
        metadata = pnd.read_csv(f"working/genomes/{meta_basename}.csv", index_col=0) 
        if all([os.path.exists(f'working/genomes/{accession}.fna') for accession in metadata['assembly_accession']]):
            metadata = metadata[metadata['local_filename'].str.endswith('assembly_stats.txt')==False]  # keep only rows for genomes
            logger.info(f"Genome assemblies already downloaded for taxids {taxids}: {len(metadata['assembly_accession'])} assemblies found. Skipping the download from NCBI.")
 
            # create metadata table and genomes dictionary: 
            get_metadata_table(logger, f'working/genomes/{meta_basename}.csv')
            create_genomes_dictionary(logger)
            response = update_metadata_manual(logger, metadata_man, source='species_to_genome')
            if response==1: return 1


            return 0    
    
          
    # download from ncbi: 
    get_genomes(logger, taxids, cores)

    
    # create the metadata table and the genomes dictionary
    get_metadata_table(logger, f'working/genomes/{meta_basename}.csv')
    create_genomes_dictionary(logger)
    response = update_metadata_manual(logger, metadata_man, source='species_to_genome')
    if response==1: return 1
    
    return 0 
    
    

def handle_manual_genomes(logger, genomes, metadata):
    
    
    # create a sub-directory without overwriting
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # create a species-to-genome dictionary
    species_to_genome = {}
    logger.debug(f"Checking the formatting of the provided -g/-genomes attribute...") 
    
    
    # check if the user specified a folder:
    if os.path.exists(genomes) and os.path.isdir(genomes):
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
    
    
    
    # Create the genomes/genomes.csv like if genomes were downloaded from NCBI.
    # Useful during plot generation.
    # Warning: the same columns are used in get_metadata_table(). But here only 2 can be filled: 'organism_name' and 'strain_isolate'.
    get_genomes_csv(source='species_to_genome')
    response = update_metadata_manual(logger, metadata, source='species_to_genome')
    if response==1: return 1
    
    
    
    return 0
    