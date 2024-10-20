import pickle


import pandas as pnd
import cobra


from ..recon import preliminary_checks
from ..recon import draft_reconstruction
from ..recon import automated_curation


from .priogapfiller import prio_gapfiller


from ..derive import derive_all


from ..commons import get_md5_string
from ..commons import get_outdir



def autopilot_command(args, logger):
        
    
    response = preliminary_checks(args, logger)
    if response == 0: return 0
    if response == 1: return 1
    
    
    response = draft_reconstruction(args, logger)
    if response == 1: return 1


    # insert the prioritized gapfiller, just before the deduplication.
    # this will OVERWRITE the panmodel in 'working/duplicates/draft_panmodel.json'
    # taht was created during draft_reconstruction()!  
    response = prio_gapfiller(logger, args.refmodel, args.refproteome, args.staining, args.mancor, args.media, args.minpanflux)
    if response == 1: return 1


    # save the panmodel md5:
    panmodel_md5 = get_md5_string('working/duplicates/draft_panmodel.json')
    with open('working/duplicates/md5.pickle', 'wb') as handle: 
        pickle.dump(panmodel_md5, handle)

        
    response = automated_curation(args, logger)
    if response == 1: return 1


    # make a SBML copy of the final draft panmodel if requested
    if args.sbml: 
        outdir = get_outdir(args.outdir)
        # load the final draft panmodel (prio-gapfilled)
        draft_panmodel = cobra.io.load_json_model(outdir + 'draft_panmodel.json')
        # create a SBML copy
        cobra.io.write_sbml_model(draft_panmodel, outdir + 'draft_panmodel.xml')


    # now we have an already gap-filled and (partially) curated panmodel.
    # apply the derivation of strain- and species-specific metabolic models:
    # load input files:
    outdir = get_outdir(args.outdir)
    panmodel = cobra.io.load_json_model( outdir + 'draft_panmodel.json')
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    report = pnd.read_csv(outdir + 'report.csv', index_col=0)
    
    if args.genbanks != '-':  # if the user provided genbanks, then there must but a 'gsnnots.csv' in output.
        # "low_memory=False" to prevent "DtypeWarning: Columns ... have mixed types."
        gannots = pnd.read_csv(outdir + 'gannots.csv', index_col=0, low_memory=False)
        gannots = gannots.astype(str)  # 'ncbigene' is still formatted as number
        gannots = gannots.replace('nan', None)
        #gannots['refseq'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in gannots['refseq'] ]
        gannots['ncbigene'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in gannots['ncbigene'] ]
    else: gannots = None

    
    response = derive_all(logger, outdir, args.cores, panmodel, pam, report, gannots, args.media, args.minflux, args.biolog, args.sbml, False, args.nofig, args.aux, args.cnps)   # False is for skipgf=False.
    if response == 1: return 1
    
    
    return 0
