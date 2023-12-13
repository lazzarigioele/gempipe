import argparse
import sys
import multiprocessing 
import logging 
from logging.handlers import QueueHandler
import traceback



from .commons import funcA, funcB

from .recon import recon_command
from .derive import derive_command


    


def main(): 
    

    # create the command line arguments:
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(title='gempipe subcommands', dest='subcommand', help='', required=True)
    
    
    # subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', help='Reconstruct a draft pan-model and a PAM.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    recon_parser.add_argument("-p", "--processes", metavar='', type=int, default=1, help="Number of parallel processes to use.")
    recon_parser.add_argument("-o", "--overwrite", action='store_true', help="Delete the working/ directory at the startup.")
    recon_parser.add_argument("-t", "--taxids", metavar='', type=str, default='-', help="Taxids of the species to model (comma separated, for example '252393,68334').")
    recon_parser.add_argument("-g", "--genomes", metavar='', type=str, default='-', help="Input genome files or folder containing the genomes (see documentation).")

    
    # subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', help='Derive strain- and species-specific models.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    derive_parser.add_argument("-p", "--processes", metavar='', type=int, help="How many parallel processes to use.")
   

    # check the inputted subcommand, automatic sys.exit(1) if a bad subprogram was specied. 
    args = parser.parse_args()
    
    
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
    logger.info('Welcome to gempipe! Launching the pipeline...')


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
    queue.put(None) # send the sentinel message
    logger_process.join() # wait for all logs to be digested
    if response == 1: sys.exit(1)
    else: sys.exit(0) # exit without errors
        
        
        
if __name__ == "__main__":
    main()