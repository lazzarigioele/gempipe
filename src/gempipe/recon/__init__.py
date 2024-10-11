import os
import shutil
import subprocess
import pickle


import pandas as pnd


import cobra


from .getgenomes import download_genomes
from .getgenomes import handle_manual_genomes
from .extractcds import extract_cds
from .extractcds import handle_manual_proteomes
from .handlegbff import handle_manual_genbanks
from .filtergenomes import filter_genomes
from .clustercds import compute_clusters
from .rec_masking import recovery_masking
from .rec_broken import recovery_broken
from .rec_overlap import recovery_overlap
from .funcannot import func_annot
from .networkrec import network_rec
from .tcdbing import tcdbing_main
from .reciprocalhits import perform_brh
from .reciprocalhits import convert_reference
from .refexpansion import ref_expansion
from .pimp import denovo_annotation
from .duplicates import solve_duplicates
from .reporting import create_panmodel_proteome
from .reporting import create_report
from .reporting import create_recon_plots


from ..commons import get_md5_string
from ..commons import get_outdir



def preliminary_checks(args, logger):
    
    
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



def draft_reconstruction(args, logger):
    
    
    # assure the output directory
    outdir = get_outdir(args.outdir)
    
    
    ### PART 1. Obtain the preoteomes. 
    
    if args.genbanks != '-':
        # handle the manually defined genbanks: 
        response = handle_manual_genbanks(logger, args.genbanks, outdir, args.metadata)
        if response == 1: return 1
    
    elif args.proteomes != '-':
        # handle the manually defined proteomes: 
        response = handle_manual_proteomes(logger, args.proteomes, args.metadata)
        if response == 1: return 1
    
    elif args.genomes != '-':
        # handle the manually defined genomes: 
        response = handle_manual_genomes(logger, args.genomes, args.metadata)
        if response == 1: return 1
    
        # extract the CDSs from the genomes:
        response = extract_cds(logger, args.cores, outdir, args.nofig)
        if response == 1: return 1    
    
        # filter the genomes based on technical/biological metrics:
        response = filter_genomes(logger, args.cores, args.buscodb, args.buscoM, args.buscoF, args.ncontigs, args.N50, outdir, args.nofig)
        if response == 1: return 1  
    
    elif args.taxids != '-':
        # download the genomes according to the specified taxids: 
        response = download_genomes(logger, args.taxids, args.cores, args.metadata)
        if response == 1: return 1 
    
        # extract the CDSs from the genomes:
        response = extract_cds(logger, args.cores, outdir, args.nofig)
        if response == 1: return 1 
    
        # filter the genomes based on technical/biological metrics:
        response = filter_genomes(logger, args.cores, args.buscodb, args.buscoM, args.buscoF, args.ncontigs, args.N50, outdir, args.nofig)
        if response == 1: return 1  
    
    else:
        logger.error("Please specify the input species taxids (-t/--taxids) or the input genomes (-g/--genomes) or the input proteomes (-p/--proteomes) or the input genbanks (-gb/--genbanks).")
        return 1
    
    
    ### PART 2. Clustering. 
    
    # cluster the aminoacid sequences according to sequence similarity. 
    response = compute_clusters(logger, args.cores)
    if response == 1: return 1 


    ### PART 3. Gene recovery.
    gene_recovery = True
    if (args.proteomes != '-' or args.genbanks != '-') and args.norec == False:  # warning if starting from proteomes
        logger.warning("gempipe gives its best when starting from genomes. Starting from proteomes/genbanks will skip the gene recovery modules.")
        gene_recovery = False
        
    if args.norec:
        gene_recovery = False
    
    if gene_recovery: 
        # Recovery 1: search for proteins broken in two
        response = recovery_broken(logger, args.cores)
        if response == 1: return 1
        
        # Recovery 2: search missing genes after masking the genome 
        response = recovery_masking(logger, args.cores)
        if response == 1: return 1 
        
        # Recovery 3: search for overlapping genes
        response = recovery_overlap(logger, args.cores)
        if response == 1: return 1
    
    # define the final PAM in the current directory: 
    if gene_recovery:
        shutil.copyfile('working/rec_overlap/pam.csv', outdir + 'pam.csv')
    else: shutil.copyfile('working/clustering/pam.csv', outdir + 'pam.csv')
    
    
    ### PART 4. Reconstruction of the reference-free reaction network.
    
    # perform functional annotation
    response = func_annot(logger, args.cores, outdir, args.dbs, args.dbmem)
    if response == 1: return 1

    # define the final annotation in the current directory: 
    annotation = pnd.read_csv('working/annotation/pan.emapper.annotations', sep='\t', comment='#', header=None)
    annotation.columns = 'query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs'.split('\t')
    annotation = annotation.set_index('query', drop=True, verify_integrity=True)
    annotation.to_csv(outdir + 'annotation.csv')

    # perform the reaction network reconstruction
    response = network_rec(logger, args.cores, args.staining, args.identity, args.coverage, args.refmodel, args.refproteome)
    if response == 1: return 1
    
    
    ### PART 5. Eventual reference-based reconstruction.
    
    if args.refmodel != '-' and args.refproteome != '-':
        
        # compute the best reciprocal hits for all the strains:
        response = perform_brh(logger, args.cores, args.refproteome)
        if response == 1: return 1
        
        # convert the reference model's genes to clusters: 
        response = convert_reference(logger, args.refmodel, args.refproteome, gene_recovery, args.refspont)
        if response == 1: return 1
    
        # expand the reference model with new reactions coming from the reference-free recon.
        response = ref_expansion(logger, args.refmodel, args.mancor, args.identity, args.coverage)
        if response == 1: return 1
    
    
    # define the final draft pan-model in the duplicates/ directory: 
    os.makedirs('working/duplicates/', exist_ok=True)  # without overwriting
    if args.refmodel != '-' and args.refproteome != '-':
        shutil.copyfile('working/expansion/draft_panmodel.json', 'working/duplicates/draft_panmodel.json')
    else: shutil.copyfile(f'working/free/draft_panmodel_{args.identity}_{args.coverage}.json', 'working/duplicates/draft_panmodel.json')
    
    
    ### PART X: Eventual TCDB-based transporters recon (experimental)
    
    if args.tcdb:  # experimental feature
        # try to reconstruct transport reactions using the tcdb. 
        response = tcdbing_main(logger, args.cores, args.staining)
        if response == 1: return 1
        
        
    return 0 



def automated_curation(args, logger):
    
    
    # assure the output directory: 
    outdir = get_outdir(args.outdir)
    
    
    # PART 6. Automated curation
    
    # perform de-novo annotation with metanetx 4.4 (aka pimp_my_model)
    response = denovo_annotation(logger)
    if response == 1: return 1


    #if args.refmodel != '-': 
    if args.dedup :
        # solve duplicate metabolites and reactions using mnx annotation
        response = solve_duplicates(logger, args.identity, args.coverage, args.refmodel, args.mancor)
        if response == 1: return 1

        # save the final panmodel
        shutil.copyfile(f'working/duplicates/draft_panmodel_da_dd.json', outdir + 'draft_panmodel.json')

    else:
        shutil.copyfile(f'working/duplicates/draft_panmodel_da.json', outdir + 'draft_panmodel.json')
    
        
    # PART 8. Create the reference panmodel proteome
    
    response = create_panmodel_proteome(logger, outdir)
    if response == 1: return 1
        
    
    # PART 7. Report
    
    # create a report table
    response = create_report(logger, outdir)
    if response == 1: return 1

    # create plots showing >preliminar< reconstruction metrics
    create_recon_plots(logger, outdir, args.cores, args.nofig)
    if response == 1: return 1

    
    return 0



def recon_command(args, logger):

    
    # overwrite, list databases, check staining ...
    response = preliminary_checks(args, logger)
    if response == 0: return 0
    if response == 1: return 1
    
    
    # produce a reference-free / reference-based reconstruction
    response = draft_reconstruction(args, logger)
    if response == 1: return 1


    # save the panmodel md5:
    panmodel_md5 = get_md5_string('working/duplicates/draft_panmodel.json')
    with open('working/duplicates/md5.pickle', 'wb') as handle: 
        pickle.dump(panmodel_md5, handle)
    
    
    # apply some automated curation (duplicates removal, ...)
    response = automated_curation(args, logger)
    if response == 1: return 1


    # make a SBML copy of the final draft panmodel if requested
    if args.sbml: 
        outdir = get_outdir(args.outdir)
        # load the final draft panmdoel
        draft_panmodel = cobra.io.load_json_model(outdir + 'draft_panmodel.json')
        # create a SBML copy
        cobra.io.write_sbml_model(draft_panmodel, outdir + 'draft_panmodel.xml')
    
    
    return 0
    
    
        
        
        
    
    