import os
import shutil
import sys
import logging
from logging.handlers import QueueHandler


from .download import get_genomes




def recon_command(args, logger):

    
    # overwrite if requested:
    if os.path.exists('working/'):
        logger.info("Found a previously created ./working/ directory.")
    if args.overwrite:
        logger.info("Ereasing the ./working/ directory as requested.")
        shutil.rmtree('working/')  
        os.makedirs('working/')
    
    
    # control of the input:
    if args.taxids == '-' or args.taxids == None or args.taxids == '': 
        logger.error("Please specify the species taxids.")
        return 1
        
        
    # download the genomes: 
    get_genomes(args.taxids)
    