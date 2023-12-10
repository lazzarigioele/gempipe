import argparse



from .moduleA import funcA
from .moduleB import funcB



def recon_command(args):
    print(args)
    
    
def derive_command(args):
    print(args)





def main(): 
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(title='', dest='subcommand', help='')

    # Subparser for the 'recon' command
    recon_parser = subparsers.add_parser('recon', help='')
    recon_parser.add_argument('arg1', help='')
    recon_parser.add_argument('arg2', help='')

    # Subparser for the 'derive' command
    derive_parser = subparsers.add_parser('derive', help='')
    derive_parser.add_argument('arg1', help='')
    derive_parser.add_argument('arg2', help='')

    
    args = parser.parse_args()
    if args.subcommand == 'recon':
        recon_command(args)
    elif args.subcommand == 'derive':
        derive_command(args)
    else:
        print("Invalid subcommand. Use 'recon' or 'derive'. EEEE")
        
        
        
if __name__ == "__main__":
    
    
    main()