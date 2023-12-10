import argparse



from .commons import funcA, funcB

from .recon import recon_command
from .derive import derive_command




def main(): 
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(title='gempipe subcommands', dest='subcommand', help='', required=True)
    
    # Subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', help='Reconstruct a draft pan-model and a PAM.')
    recon_parser.add_argument("-p", "--processes", metavar='', type=int, help="How many parallel processes to use.")
    recon_parser.add_argument("-a", "--optionA", metavar='', help="Option A for recon")
    recon_parser.add_argument("-b", "--optionB", metavar='', help="Option B for recon")
    recon_parser.add_argument("-c", "--optionC", metavar='', help="Option C for recon")

    # Subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', help='Derive strain- and species-specific models.')
    # optional
    derive_parser.add_argument("-p", "--processes", metavar='', type=int, help="How many parallel processes to use.")
    derive_parser.add_argument("-x", "--optionX", metavar='', help="Option X for derive")
    derive_parser.add_argument("-y", "--optionY", metavar='', help="Option Y for derive")
    # positional
    derive_parser.add_argument("N", help="Option N for derive")
   

    # Call the appropriate function based on the sub-command
    args = parser.parse_args()
    if args.subcommand == 'recon':
        recon_command(args)
    elif args.subcommand == 'derive':
        derive_command(args)
    else:
        print("Invalid subcommand.")
        
        
        
if __name__ == "__main__":

    main()