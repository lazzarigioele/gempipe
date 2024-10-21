import os


import pandas as pnd
import cobra 


from .strain import derive_strain_specific
from .filler import strain_species_filler
from .biolog import strain_biolog_tests
from .screening import strain_auxotrophies_tests
from .screening import strain_cnps_tests
from .species import derive_rpam
from .species import derive_species_specific
from .reporting import create_derive_plots



def derive_all(logger, outdir, cores, panmodel, pam, report, gannots, media_filepath, minflux, biolog, sbml, skipgf, nofig, aux, cnps):
    
    
    ### PART 1: derive strain-specific models
    response = derive_strain_specific(logger, outdir, cores, panmodel, pam, report, gannots, sbml)
    if response == 1: return 1


    ### PART 2: gap-fill strain-specific models
    if skipgf: logger.info("Skipping the strain-specific gap-filling step (--skipgf)...")
    if not skipgf:
        response = strain_species_filler(logger, outdir, cores, panmodel, media_filepath, minflux, sbml, level='strain')
        if response == 1: return 1
    

    ### PART 2bis: simulations
    if biolog:
        response = strain_biolog_tests(logger, outdir, cores, pam, panmodel, skipgf)
        if response == 1: return 1
    
    if aux:
        response = strain_auxotrophies_tests(logger, outdir, cores, pam, skipgf)
        if response == 1: return 1
    
    if cnps:
        response = strain_cnps_tests(logger, outdir, cores, pam, panmodel, skipgf)
        if response == 1: return 1
    
     
    ### PART 3: derive species-specific models
    response = derive_rpam(logger, outdir, cores, panmodel, skipgf)
    if response == 1: return 1
    
    response = derive_species_specific(logger, outdir, cores, panmodel, sbml)
    if response == 1: return 1

    if skipgf: logger.info("Skipping the species-specific gap-filling step (--skipgf)...")
    if not skipgf:
        response = strain_species_filler(logger, outdir, cores, panmodel, media_filepath, minflux, sbml, level='species')
        if response == 1: return 1
    
    
    ### PART 4: make some plots
    create_derive_plots(logger, outdir, nofig)
    if response == 1: return 1
    

    return 0



def derive_command(args, logger):
    
    
    # check the existence of the input files:
    if args.inpanmodel == '-' or args.inpam == '-' or args.inreport == '-':
        logger.error("Please specify the input pan-model (-im/--inpanmodel), PAM (-ip/--inpam) and report (-ir/--inreport).")
        return 1
    else:  # all 3 parameters were set
        if not os.path.exists(args.inpanmodel):
            logger.error(f"The specified path for input pan-model (-im/--inpanmodel) does not exist: {args.inpanmodel}.")
            return 1
        if not os.path.exists(args.inpam):
            logger.error(f"The specified path for input PAM (-ip/--inpam) does not exist: {args.inpam}.")
            return 1
        if not os.path.exists(args.inreport):
            logger.error(f"The specified path for input report (-ir/--inreport) does not exist: {args.inreport}.")
            return 1
        
        
    if args.ingannots != '-':
        if not os.path.exists(args.ingannots):
            logger.error(f"The specified path for input genes annotation table (-ig/--ingannots) does not exist: {args.ingannots}.")
            return 1
    
    
    # load input files
    logger.info("Loading input files...")
    panmodel = cobra.io.load_json_model(args.inpanmodel)
    pam = pnd.read_csv(args.inpam, index_col=0)
    report = pnd.read_csv(args.inreport, index_col=0)
    
    if args.ingannots != '-':
        # "low_memory=False" to prevent "DtypeWarning: Columns ... have mixed types."
        gannots = pnd.read_csv(args.ingannots, index_col=0, low_memory=False)
        gannots = gannots.astype(str)  # 'ncbigene' is still formatted as number
        gannots = gannots.replace('nan', None)
        #gannots['refseq'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in gannots['refseq'] ]
        gannots['ncbigene'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in gannots['ncbigene'] ]
    else: gannots = None
    
    
    # create the main output directory: 
    outdir = args.outdir
    if outdir.endswith('/') == False: outdir = outdir + '/'
    os.makedirs(outdir, exist_ok=True)
    
    
    logger.info("Deriving strain- and species-specific metabolic models...")
    response = derive_all(logger, outdir, args.cores, panmodel, pam, report, gannots, args.media, args.minflux, args.biolog, args.sbml, args.skipgf, args.nofig, args.aux, args.cnps)
    if response == 1: return 1
    
    
    return 0
