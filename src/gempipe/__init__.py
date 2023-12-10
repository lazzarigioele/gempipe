import argparse



from .commons import funcA, funcB

from .recon import recon_command
from .derive import derive_command




def main(): 
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(title='gempipe subcommands', dest='subcommand', help='')

    # Subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', help='Reconstruct a draft pan-model and a PAM.')
    recon_parser.add_argument('--arg1', help='Description for --arg1 argument.')
    recon_parser.add_argument('--arg2', help='Description for --arg2 argument.')

    # Subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', help='Derive strain- and species-specific models.')
    derive_parser.add_argument('--arg1', help='Description for --arg1 argument.')
    derive_parser.add_argument('--arg2', help='Description for --arg2 argument.')

    
    args = parser.parse_args()
    if args.subcommand == 'recon':
        recon_command(args)
    elif args.subcommand == 'derive':
        derive_command(args)
    else:
        print("Invalid subcommand.")
        
        
        
if __name__ == "__main__":
    
    
    main()