import os
import shutil



from .download import get_genomes
from .download import get_metadata_table
from .download import handle_manual_genomes


from .annotatecds import extract_cds
from .annotatecds import handle_manual_proteomes



def recon_command(args, logger):

    
    # overwrite if requested:
    if os.path.exists('working/'):
        logger.info("Found a previously created ./working/ directory.")
    if args.overwrite:
        logger.info("Ereasing the ./working/ directory as requested (-o/--overwrite).")
        shutil.rmtree('working/')  
        os.makedirs('working/')
        
    
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
    
    else:
        logger.error("Please specify the species taxids (-t/--taxids) or the input genomes (-g/--genomes) or the input proteomes (-p/--proteomes).")
        return 1
    