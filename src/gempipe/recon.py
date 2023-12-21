import os
import shutil
import subprocess



from .getgenomes import download_genomes
from .getgenomes import handle_manual_genomes
from .extractcds import extract_cds
from .extractcds import handle_manual_proteomes
from .filtergenomes import filter_genomes
from .clustercds import compute_clusters
from .rec_masking import recovery_masking
from .rec_broken import recovery_broken
from .rec_overlap import recovery_overlap



def recon_command(args, logger):

    
    # overwrite if requested:
    if os.path.exists('working/'):
        logger.info("Found a previously created ./working/ directory.")
        if args.overwrite:
            logger.info("Ereasing the ./working/ directory as requested (--overwrite).")
            shutil.rmtree('working/')  
    os.makedirs('working/', exist_ok=True)
    os.makedirs('working/logs/', exist_ok=True)
        
        
    # check if the user required the list of databases: 
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
        response = download_genomes(logger, args.taxids, args.cores)
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
    
    
    
    ### PART 2. Clustering. 
    
    # cluster the aminoacid sequences according to sequence similarity. 
    response = compute_clusters(logger, args.cores)
    if response == 1: return 1 



    ### PART 3. Gene recovery.
    
    if args.genomes != '-' or args.taxids != '-': 
        # Recovery 1: search for proteins broken in two
        response = recovery_broken(logger, args.cores)
        if response == 1: return 1
        
        # Recovery 2: search missing genes after masking the genome 
        response = recovery_masking(logger, args.cores)
        if response == 1: return 1 
        
        # Recovery 3: search for overlapping genes
        response = recovery_overlap(logger, args.cores)
        if response == 1: return 1
        
        
        
    # warning if starting from proteomes
    if args.proteomes != '-': 
        logger.warning("gempipe gives its best when starting from genomes. Starting from proteomes will skip the gene recovery modules.")
    