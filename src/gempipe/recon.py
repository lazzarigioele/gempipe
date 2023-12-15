import os
import shutil
import subprocess



from .download import get_genomes
from .download import get_metadata_table
from .download import handle_manual_genomes


from .annotatecds import extract_cds
from .annotatecds import handle_manual_proteomes


from .filtergenomes import filter_genomes



def recon_command(args, logger):

    
    # overwrite if requested:
    if os.path.exists('working/'):
        logger.info("Found a previously created ./working/ directory.")
        if args.overwrite:
            logger.info("Ereasing the ./working/ directory as requested (-o/--overwrite).")
            shutil.rmtree('working/')  
    os.makedirs('working/', exist_ok=True)
    os.makedirs('working/logs/', exist_ok=True)
        
        
    # check is the user required the list of databases: 
    if args.buscodb == 'show': 
        logger.info("Creating the temporary ./busco_downloads/ directory...")
        command = f"""busco --list-datasets"""
        process = subprocess.Popen(command, shell=True)
        process.wait()
        shutil.rmtree('busco_downloads/') 
        logger.info("Deleted the temporary ./busco_downloads/ directory.")
        return 0
        
    
    # check inputted gram staining 
    if args.staining != 'pos' and args.staining != 'neg': 
        logger.error("Gram staining (-s/--staining) must be either 'pos' or 'neg'.")
        return 1
    
    
    
    ### PART 1. Obtain the preoteomes. 
    
    if args.proteomes != '-':
        # handle the manually defined proteomes: 
        response = handle_manual_proteomes(logger, args.proteomes)
        if response == 1: return 1
    
    elif args.genomes != '-':
        # handle the manually defined genomes: 
        response = handle_manual_genomes(logger, args.genomes)
        if response == 1: return 1
    
        # extract the CDSs from the genomes:
        response = extract_cds(logger, args.cores)
        if response == 1: return 1    
    
        # filter the genomes based on technical/biological metrics:
        response = filter_genomes(logger, args.cores, args.buscodb, args.buscoM, args.ncontigs, args.N50)
        if response == 1: return 1  
    
    elif args.taxids != '-':
        # download the genomes according to the specified taxids: 
        response = get_genomes(logger, args.taxids, args.cores)
        if response == 1: return 1
    
        # get the metadata table:
        response = get_metadata_table(logger)
        if response == 1: return 1
    
        # extract the CDSs from the genomes:
        response = extract_cds(logger, args.cores)
        if response == 1: return 1 
    
        # filter the genomes based on technical/biological metrics:
        response = filter_genomes(logger, args.cores, args.buscodb, args.buscoM, args.ncontigs, args.N50)
        if response == 1: return 1  
    
    else:
        logger.error("Please specify the species taxids (-t/--taxids) or the input genomes (-g/--genomes) or the input proteomes (-p/--proteomes).")
        return 1
    