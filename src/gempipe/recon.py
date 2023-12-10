import os
import shutil


from .download import get_genomes



def recon_command(args):

    
    # overwrite if requested:
    if args.overwrite:
        if os.path.exists('working/'): shutil.rmtree('working/')  
    os.makedirs('working/', exist_ok=True)
    
    
    # control of the input:
    if args.taxids == '-' or args.taxids == None or args.taxids == '': 
        print("Please specify the species taxids.")
        
        
    # download the genomes: 
    get_genomes(args.taxids)
    