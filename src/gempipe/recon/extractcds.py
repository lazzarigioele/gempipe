import pickle
import os
import subprocess
import multiprocessing
import itertools
import shutil
import glob


import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


from ..commons import chunkize_items
from ..commons import load_the_worker
from ..commons import gather_results

from ..commons import get_genomes_csv
from ..commons import update_metadata_manual




            
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
        
        
    # convert gff files to dataframe: 
    coords_df = []
    with open(f'working/proteomes/{accession}.gff', 'r') as r_handler: 
        file = r_handler.read()
        file = file.split('##FASTA', 1)[0] # not interested in contigs
        for line in file.split('\n'):
            if line.startswith('#'): continue
            if line=='': continue
            row = line.split('\t')
            row_dict = {key: value for key, value in zip(['contig', 'col2', 'col3', 'start', 'end', 'col6', 'strand', 'col8', 'attributes'], row)}
            # extract attributes from col9:
            for attribute in row_dict['attributes'].split(';'):
                attribute_id, attribute_value = attribute.split('=', 1)
                row_dict[attribute_id] = attribute_value
            # add the accession column:
            row_dict['accession'] = accession
            coords_df.append(row_dict)
    coords_df = pnd.DataFrame.from_records(coords_df)
    
    
    # format the dataframe and save
    coords_df = coords_df[['ID', 'accession', 'contig', 'strand', 'start', 'end']]
    coords_df.to_csv(f'working/coordinates/{accession}.csv')
        
        
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
    # gff files could be useful for pangenome analysis like Raory/Panaroo:
    shutil.move(f'working/proteomes/{accession}.gff', f'working/gff/{accession}.gff')
        
    
    # return a row for the dataframe
    return [{'accession': accession, 'completed': True}]



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
    
    
    
def create_seq_to_coords(logger):
    
    
    # load the previously created species_to_proteome: 
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)

    
    # create a list of accessions to parse: 
    accessions = set()
    for species in species_to_proteome.keys(): 
        for proteome in species_to_proteome[species]:
            basename = os.path.basename(proteome)
            accession, _ = os.path.splitext(basename)
            accessions.add(accession)
            
    
    # quick check of the available dict to save time:
    if os.path.exists(f'working/coordinates/seq_to_coords.pickle'):
        accessions_available = set()
        with open('working/coordinates/seq_to_coords.pickle', 'rb') as handler:
            seq_to_coords = pickle.load(handler)
        for attribs in seq_to_coords.values():
            accessions_available.add(attribs['accession'])
        if accessions == accessions_available:
            logger.debug(f"Found a good sequence-to-coordinates dictionary file already available.")
            return
    
        
    # create the dict seq-to-coords
    seq_to_coords = {}
    for accession in accessions: 
        coords_df = pnd.read_csv(f'working/coordinates/{accession}.csv')
        for index, row in coords_df.iterrows(): 
            if type(row['ID']) != str:
                continue
            # 'str' needed as some accessions could be just numbers
            seq_to_coords[row['ID']] = {'accession': str(row['accession']), 'contig': str(row['contig']), 'strand': str(row['strand']), 'start': row['start'], 'end': row['end']}
    
    
    # save the dictionary to disk: 
    with open('working/coordinates/seq_to_coords.pickle', 'wb') as file:
        pickle.dump(seq_to_coords, file)
    logger.debug(f"Saved the sequence-to-coordinates dictionary to file: ./working/coordinates/seq_to_coords.pickle.")
    


def figure_cds(logger, outdir):
    
    
    logger.info("Producing figure for extracted CDSs in {outdir}/figures/n_cds.png...")
    
    # create summary dataframe  (should be equivalent to 'all_df_combined'): 
    prodigal_summary = []
    for file in glob.glob(f"working/proteomes/*.faa"):
        accession = file.replace('working/proteomes/', '').replace('.faa', '')
        n_cds = open(file).read().count('>')
        prodigal_summary.append({'assembly_accession': accession, 'n_cds': n_cds})
    prodigal_summary = pnd.DataFrame.from_records(prodigal_summary)
    prodigal_summary = prodigal_summary.set_index('assembly_accession', drop=True, verify_integrity=True)
    
    # load the genomes_df to have the 'strain_isolate' and 'organism_name' columns:
    genomes_df = pnd.read_csv('working/genomes/genomes.csv', index_col=0)
    genomes_df = genomes_df.set_index('assembly_accession', drop=True, verify_integrity=True)
    
    # concat the dataframes:
    df = pnd.concat([genomes_df, prodigal_summary], axis=1)
    
    # define colors:
    df = df.set_index('strain_isolate', drop=False)
    colors = df['organism_name'].map({species: f'C{number}' for number, species in enumerate(df['organism_name'].unique())}).to_dict()
        
    # draw bars:
    fig, ax = plt.subplots()
    _ = sb.barplot(df, x='strain_isolate', y='n_cds', palette=colors, hue='strain_isolate', legend=False, ax=ax)
    
    # set tick labels
    ax.tick_params(axis='x', labelrotation=90)
    [label.set_color(colors[label.get_text()]) for label in ax.get_xticklabels()]
    
    # set legend:
    plt.legend(handles=[Patch(color=f'C{number}', label=species) for number, species in enumerate(df['organism_name'].unique())], title='', loc='center left', bbox_to_anchor=(1.05, 0.5))
    
    ax.figure.set_size_inches(0.2*len(df), 4)
    sb.despine()

    
    if len(df) <= 100:
        plt.savefig(outdir + 'figures/n_cds.png', dpi=300, bbox_inches='tight')
    else:
        logger.info("Number of genomes is >100: producing the SVG version instead {outdir}/figures/n_cds.svg...")
        plt.savefig(outdir + 'figures/n_cds.svg', bbox_inches='tight')
        
    
    
    
def extract_cds(logger, cores, outdir, nofig):
    
    
    # create sub-directory without overwriting:
    logger.info("Extracting the CDSs from the genomes...")
    os.makedirs('working/proteomes/', exist_ok=True)
    os.makedirs('working/coordinates/', exist_ok=True)
    os.makedirs('working/gff/', exist_ok=True)
    os.makedirs(outdir + 'figures/', exist_ok=True)


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
        already_computed.append(os.path.exists(f'working/coordinates/{accession}.csv'))
    if all(already_computed):
        logger.info("Found all the needed files already computed. Skipping this step.")
        # save the species_to_proteome and seq_to_coords dicts
        create_species_to_proteome(logger)
        create_seq_to_coords(logger)
        
        if not nofig:
            figure_cds(logger, outdir)
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
    all_df_combined = gather_results(results)  # all_df_combined can be ignored.
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join()  
    
    
    # save the species_to_proteome and seq_to_coords dicts
    create_species_to_proteome(logger)
    create_seq_to_coords(logger)
    
    
    if not nofig:
        figure_cds(logger, outdir)
    return 0



def handle_manual_proteomes(logger, proteomes, metadata):
    
    
    # create a species-to-genome dictionary
    species_to_proteome = {}
    logger.debug(f"Checking the formatting of the provided -p/--proteomes attribute...") 
    
    
    # check if the user specified a folder:
    if os.path.exists(proteomes) and os.path.isdir(proteomes):
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
            basename = os.path.basename(file)
            shutil.copyfile(file, 'working/proteomes/' + basename)
            copied_files.append('working/proteomes/' + basename)
        species_to_proteome[species] = copied_files
    logger.debug(f"Input proteomes copied to ./working/proteomes/.")
    logger.debug(f"Created the species-to-proteome dictionary: {str(species_to_proteome)}.") 
    
    
    # save the dictionary to disk: 
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome, file)
    logger.debug(f"Saved the species-to-proteome dictionary to file: ./working/proteomes/species_to_proteome.pickle.")
    
    
    # Create the genomes/genomes.csv like if genomes were downloaded from NCBI.
    # Useful during plot generation.
    # Warning: the same columns are used in get_metadata_table(). But here only 2 can be filled: 'organism_name' and 'strain_isolate'.
    get_genomes_csv(source='species_to_proteome')
    response = update_metadata_manual(logger, metadata, source='species_to_proteome')
    if response==1: return 1
    
    
    
    return 0



