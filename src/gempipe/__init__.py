import argparse
import sys
import multiprocessing 
import logging 
from logging.handlers import QueueHandler



from .commons import funcA, funcB

from .recon import recon_command
from .derive import derive_command


    


def main(): 
    

    # Create the command line arguments:
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(title='gempipe subcommands', dest='subcommand', help='', required=True)
    
    
    # Subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', help='Reconstruct a draft pan-model and a PAM.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    recon_parser.add_argument("-p", "--processes", metavar='', type=int, default=1, help="Number of parallel processes to use.")
    recon_parser.add_argument("-o", "--overwrite", action='store_true', help="Delete the working/ directory as first step.")
    recon_parser.add_argument("-t", "--taxids", metavar='', type=str, default='-', help="Taxids of the species to model (comma separated, for example '252393,68334').")
    recon_parser.add_argument("-g", "--genomes", metavar='', type=str, default='-', help="Path to the input genomes (comma separated, for example 'mydir/g1.fna,mydir/g2.fna,mydir/g3.fna' or 'SpA:mydir/g1.fna;SpB:mydir/g2.fna,mydir/g3.fna').")
    recon_parser.add_argument("-a", "--optionA", metavar='', help="Option A for recon")
    recon_parser.add_argument("-b", "--optionB", metavar='', help="Option B for recon")
    recon_parser.add_argument("-c", "--optionC", metavar='', help="Option C for recon")

    
    # Subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', help='Derive strain- and species-specific models.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # optional
    derive_parser.add_argument("-p", "--processes", metavar='', type=int, help="How many parallel processes to use.")
    derive_parser.add_argument("-x", "--optionX", metavar='', help="Option X for derive")
    derive_parser.add_argument("-y", "--optionY", metavar='', help="Option Y for derive")
    # positional
    derive_parser.add_argument("N", help="Option N for derive")
   

    # Check the inputted subcommand, sys.exit(1) if a bad subprogram was specied. 
    args = parser.parse_args()
    
    
    # Create a logging queue in a dedicated process.
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
    
    
    # Connect the logger for this (main) function: 
    logger = logging.getLogger('gempipe')
    logger.addHandler(QueueHandler(queue))
    logger.setLevel(logging.DEBUG) # debug (lvl 10) and up
    
    
    # Show a welcome message:
    logger.info('Welcome to gempipe! Launching the pipeline...')

    
    # choose which subcommand to lauch: 
    if args.subcommand == 'recon':
        response = recon_command(args, logger)
    if args.subcommand == 'derive':
        response = derive_command(args)


    # terminate the program
    queue.put(None) # send the sentinel message
    logger_process.join() # wait for all logs to be digested
    if response == 1: sys.exit(1)
    else: sys.exit(0) # exit without errors
        
        
        
if __name__ == "__main__":
    main()