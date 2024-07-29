import pickle
import os


import pandas as pnd
from Bio import SeqIO, SeqRecord


from ..commons import read_refmodel



def create_panmodel_proteome(logger, outdir):
    
    
    logger.info("Creating the reference proteome for the draft pan-model...")
    
    
    
    # A) from the panmodel
    
    # read the final draft panmodel
    draft_panmodel = read_refmodel(outdir + 'draft_panmodel.json')
    genes_to_report = set()
    for g in draft_panmodel.genes:
        if g.id == 'spontaneous':
            continue
        genes_to_report.add(g.id)
        
    
    # collect the reference sequences
    sr_list = []
    added = set()
    for record in SeqIO.parse('working/clustering/representatives.ren.faa', "fasta"):
        cluster, cds, accession = record.description.split(' ')
        if cluster in genes_to_report:
            sr = SeqRecord.SeqRecord(record.seq, id=cluster, description=f'{cds} {accession}')
            sr_list.append(sr)
            added.add(cluster)
            genes_to_report = genes_to_report - added
            
            
    # if all the sequences were recovered, write the fasta:
    if genes_to_report == set():
        with open(outdir + 'draft_panmodel.faa', 'w') as w_handler:
            count = SeqIO.write(sr_list, w_handler, "fasta")
        logger.debug(f"{len(added)} reference sequences written to " + outdir + 'draft_panmodel.faa' + '.')
        
        
        
    # B) from the PAM
    """
    # read the final draft panmodel
    pam = pnd.read_csv(outdir + 'pam.csv', index_col=0)
    genes_to_report = set()
    for cluster, row in pam.iterrows():
        genes_to_report.add(cluster)
        
    
    # collect the reference sequences
    sr_list = []
    added = set()
    for record in SeqIO.parse('working/clustering/representatives.ren.faa', "fasta"):
        cluster, cds, accession = record.description.split(' ')
        if cluster in genes_to_report:
            sr = SeqRecord.SeqRecord(record.seq, id=cluster, description=f'{cds} {accession}')
            sr_list.append(sr)
            added.add(cluster)
            genes_to_report = genes_to_report - added
            
            
    # if all the sequences were recovered, write the fasta:
    if genes_to_report == set():
        with open(outdir + 'draft_panproteome.faa', 'w') as w_handler:
            count = SeqIO.write(sr_list, w_handler, "fasta")
        logger.debug(f"{len(added)} reference sequences written to " + outdir + 'draft_panproteome.faa' + '.')
    """
            


def create_report(logger, outdir):
    
    
    report = []  # list of dicts, future dataframe
    
    
    # get the retained genomes/proteomes (post filtering):
    with open('working/proteomes/species_to_proteome.pickle', 'rb') as handler:
        species_to_proteome = pickle.load(handler)
        for species in species_to_proteome.keys(): 
            for proteome in species_to_proteome[species]:
                basename = os.path.basename(proteome)
                accession, _ = os.path.splitext(basename)
                
                
                # populate the table: 
                report.append({'species': species, 'accession': accession})
                
                
    # save to file
    report = pnd.DataFrame.from_records(report)
    report.to_csv(outdir + 'report.csv')
    
    
    
    return 0 
    
    