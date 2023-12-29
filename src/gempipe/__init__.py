import argparse
import sys
import multiprocessing 
import logging 
from logging.handlers import QueueHandler
import traceback
import importlib.metadata



from .recon import recon_command
from .derive import derive_command




from .flowchart import Flowchart



def main(): 
    
    
    # define the header of main- and sub-commands. 
    pub_details = 'TODO'
    header = f'gempipe v{importlib.metadata.metadata("gempipe")["Version"]}, please cite "{pub_details}".'
    
    
    # create the command line arguments:
    parser = argparse.ArgumentParser(description=header, add_help=False)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument("-v", "--version", action="version", version=f"v{importlib.metadata.metadata('gempipe')['Version']}", help="Show version number and exit.")
    subparsers = parser.add_subparsers(title='gempipe subcommands', dest='subcommand', help='', required=True)
    
    
    # subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', description=header, help='Reconstruct a draft pan-model and a PAM.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    recon_parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    recon_parser.add_argument("-v", "--version", action="version", version=f"v{importlib.metadata.metadata('gempipe')['Version']}", help="Show version number and exit.")
    recon_parser.add_argument("-c", "--cores", metavar='', type=int, default=1, help="Number of parallel processes to use.")
    recon_parser.add_argument("--overwrite", action='store_true', help="Delete the working/ directory at the startup.")
    recon_parser.add_argument("-t", "--taxids", metavar='', type=str, default='-', help="Taxids of the species to model (comma separated, for example '252393,68334').")
    recon_parser.add_argument("-g", "--genomes", metavar='', type=str, default='-', help="Input genome files or folder containing the genomes (see documentation).")
    recon_parser.add_argument("-p", "--proteomes", metavar='', type=str, default='-', help="Input proteome files or folder containing the proteomes (see documentation).")
    recon_parser.add_argument("-s", "--staining", metavar='', type=str, default='neg', help="Gram staining, 'pos' or 'neg'.")
    recon_parser.add_argument("-b", "--buscodb", metavar='', type=str, default='bacteria_odb10', help="Busco database to use ('show' to see the list of available databases).")
    recon_parser.add_argument("--buscoM", metavar='', type=str, default='2%', help="Maximum number of missing Busco's single copy orthologs (absolute or percentage).")
    recon_parser.add_argument("--ncontigs", metavar='', type=int, default=200, help="Maximum number of contigs allowed per genome.")
    recon_parser.add_argument("--N50", metavar='', type=int, default=50000, help="Minimum N50 allowed per genome.")
    recon_parser.add_argument("--identity", metavar='', type=int, default=30, help="Minimum percentage amino acidic sequence identity to use when aligning against the BiGG gene database.")
    recon_parser.add_argument("--coverage", metavar='', type=int, default=70, help="Minimum percentage coverage to use when aligning against the BiGG gene database.")
    recon_parser.add_argument("-rm", "--ref_model", metavar='', type=str, default='-', help="Model to be used as reference.")
    recon_parser.add_argument("-rp", "--ref_proteome", metavar='', type=str, default='-', help="Proteome to be used as reference.")
    recon_parser.add_argument("-o", "--outdir", metavar='', type=str, default='./', help="Main output directory (will be created if not existing).")
    

    
    # subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', description=header, help='Derive strain- and species-specific models.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    derive_parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    derive_parser.add_argument("-v", "--version", action="version", version=f"v{importlib.metadata.metadata('gempipe')['Version']}", help="Show version number and exit.")
    derive_parser.add_argument("-c", "--cores", metavar='', type=int, help="How many parallel processes to use.")
    derive_parser.add_argument("-o", "--outdir", metavar='', type=str, default='./', help="Main output directory (will be created if not existing).")
    

    # check the inputted subcommand, automatic sys.exit(1) if a bad subprogram was specied. 
    args = parser.parse_args()
    
    
    # set the multiprocessing context
    multiprocessing.set_start_method('fork') 
    
    
    # create a logging queue in a dedicated process.
    def logger_process_target(queue):
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger = logging.getLogger('gempipe')
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG) # debug (lvl 10) and up
        while True:
            message = queue.get() # block until a new message arrives
            if message is None: # sentinel message to exit the loop
                break
            logger.handle(message)
    queue = multiprocessing.Queue()
    logger_process = multiprocessing.Process(target=logger_process_target, args=(queue,))
    logger_process.start()
    
    
    # connect the logger for this (main) process: 
    logger = logging.getLogger('gempipe')
    logger.addHandler(QueueHandler(queue))
    logger.setLevel(logging.DEBUG) # debug (lvl 10) and up
    
    
    # show a welcome message:
    print('\n' + header + '\n', file=sys.stdout)
    logger.info("Welcome to gempipe! Launching the pipeline...")


    try: 
        # choose which subcommand to lauch: 
        if args.subcommand == 'recon':
            response = recon_command(args, logger)
        if args.subcommand == 'derive':
            response = derive_command(args)
    except: 
        # show the error stack trace for this un-handled error: 
        response = 1
        logger.error(traceback.format_exc())


    # Terminate the program:
    if response == 1: 
        queue.put(None) # send the sentinel message
        logger_process.join() # wait for all logs to be digested
        sys.exit(1)
    else: 
        # show a bye message
        logger.info("gempipe terminated without errors!")
        queue.put(None) # send the sentinel message
        logger_process.join() # wait for all logs to be digested
        print('\n' + header + '\n', file=sys.stdout)
        sys.exit(0) # exit without errors
        
        
        
if __name__ == "__main__":
    main()
    
    # some command line examples to test the pipeline
    # gempipe recon -s neg -t 252393,68334 -b enterobacterales_odb10 --ncontigs 2000 --N50 5000 --buscoM 10% -c 8
    # gempipe recon -s pos -g testing/plantarum_Siezen2010 -b lactobacillales_odb10 -c 8 -rm testing/from_mendoza/iLP728.xml -rp testing/from_mendoza/protein_fasta.faa
    # gempipe recon -s pos -g testing/plantarum_Siezen2010 -b lactobacillales_odb10 -c 8 -rm testing/from_mendoza/iLP728_improved.json -rp testing/from_mendoza/protein_fasta.faa