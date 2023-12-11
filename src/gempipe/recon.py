import os
import shutil



from .download import get_genomes
from .download import get_metadata_table
from .download import handle_manual_genomes




def recon_command(args, logger):

    
    # overwrite if requested:
    if os.path.exists('working/'):
        logger.info("Found a previously created ./working/ directory.")
    if args.overwrite:
        logger.info("Ereasing the ./working/ directory as requested (-o/--overwrite).")
        shutil.rmtree('working/')  
        os.makedirs('working/')
        
        

    
    if args.genomes != '-':
        # handle the manually defined genomes: 
        response = handle_manual_genomes(logger, args.genomes)
        if response == 1: return 1
    
    elif args.taxids != '-':
        # download the genomes according to the specified taxids: 
        response = get_genomes(logger, args.taxids, args.processes)
        if response == 1: return 1
    
        # get the metadata table:
        response = get_metadata_table(logger)
        if response == 1: return 1
    
    else:
        logger.error("Please specify the species taxids (-t/--taxids) or the input genomes (-g/--genomes).")
        return 1
    