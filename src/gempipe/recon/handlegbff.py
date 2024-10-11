import os
import glob
import subprocess
import pickle


import pandas as pnd
from Bio import SeqIO, SeqRecord, Seq


from ..commons import get_genomes_csv
from ..commons import update_metadata_manual





def get_sequences_and_source(gbff_file, ):

    # use SeqIO.read() if you expect only one record in the GenBank file.
    # alternatively, use SeqIO.parse() if there are multiple records in the file.
    records = SeqIO.parse(gbff_file, 'genbank')
                              
    sequences = []
    sourcedb = None

    for record in records:

        sequence_id = record.id
        sequence_description = record.description
        sequence_seqlen = len(record.seq)
        
        
        # append sequence id
        sequences.append(sequence_id)
        
        
    # take the first and search for an underscore (_). 
    # if found, we assume the gbff is from refseq, not genbank.
    source_db = 'refseq' if '_' in sequences[0] else 'genbank'
                              
    return sequences, source_db



def get_uniprot_table(accession, sequences, source_db, rformat):
    
    
    # choose the crossref db (embl can be seen as a replacement for genbank)
    if source_db=='genbank': 
        crossref = 'embl'
    else: crossref = 'refseq'
    
    
    # adjust sequences if from genbank: 
    if crossref=='embl':
        sequences = [i.rsplit('.', 1)[0] for i in sequences]


    # define the query string based on the format:
    query_string = '%28' + '+OR+'.join([f'%28xref%3A{crossref}-{i}%29' for i in sequences]) + '%29'

    if rformat == 'json':  
        # 'json' parsing still to implement. 'tsv' preferred because it weights much less.
        query_string = f'https://rest.uniprot.org/uniprotkb/stream?format={rformat}&query={query_string}'

    if rformat == 'tsv':  
        query_fields = ['accession', 'gene_oln', 'xref_refseq', 'xref_embl', 'xref_kegg', 'xref_geneid',]
        query_fields = '%2C'.join(query_fields) 
        query_string = f'https://rest.uniprot.org/uniprotkb/stream?fields={query_fields}&format={rformat}&query={query_string}'


    # define the tmp output text file:  
    outpath = f'working/gannots/response_{accession}.tsv'

    
    # perform the query: 
    command = f'''wget -O {outpath} --quiet "{query_string}"'''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()


    response = pnd.read_csv(outpath, sep='\t')
    os.remove(outpath)  # clean up, table is in ram.
    return response
    
    
    
def get_features(accession, gbff_file):
    
    
    # use SeqIO.read() if you expect only one record in the GenBank file.
    # alternatively, use SeqIO.parse() if there are multiple records in the file.
    records = SeqIO.parse(gbff_file, 'genbank')
                              

    results_df = []
    for record in records:

        sequence_id = record.id
        sequence_description = record.description
        sequence_seqlen = len(record.seq)
        
        
        # access features:
        for feature in record.features:
            if feature.type == 'CDS':
                try: location = feature.location
                except: location = None
                try: gene_symbol = feature.qualifiers['gene'][0]
                except: gene_symbol = None
                try: aa_seq = feature.qualifiers['translation'][0]
                except: 
                    try: aa_seq = str(feature.extract(record.seq).translate())  # reverse complement automatically handled
                    except: aa_seq = None
                try: locus_tag = feature.qualifiers['locus_tag'][0]
                except: locus_tag = None
                try: old_locus_tag = feature.qualifiers['old_locus_tag'][0]
                except: old_locus_tag = None
                if locus_tag != None and (type(locus_tag)==int or str(locus_tag).isdigit()):  
                    # gene nume is just a number (the genbank may come from prodigal/bactabolize like pipelines)
                    # probably the gene ID will conflict with genes from other strains. 
                    # Therefore we prepend its accession:
                    locus_tag = accession + '_' + record.id + '_' + str(locus_tag)
                    locus_tag = locus_tag.replace('.', '').replace('-', '')  # remove undesired chars
                    
                # Note on the 'protein_id' field: 
                # if the .gbff is from GenBank, then 'protein_id' is the actual protein (example: CCC77580.1)
                # if the .gbff is from RefSeq, then 'protein_id' is a NR protein (prefix WP_)
                try: protein_id = feature.qualifiers['protein_id'][0]
                except: protein_id = None
                
                results_df.append({
                    'accession': accession, 'seqid': record.id, 'gene_symbol': gene_symbol, # 'location': location,  
                    'locus_tag': locus_tag, 'old_locus_tag': old_locus_tag, 'protein_id': protein_id})
    
    
    results_df = pnd.DataFrame.from_records(results_df)
    return results_df
    
    
    
def get_gannots_tabular(accession, gbff_file, outfolder):
    
    
    # get accession list and source database (genbank/refseq):
    sequences, source_db = get_sequences_and_source(gbff_file, )
    
    # get the uniprot crossref table:
    up_table = get_uniprot_table(accession, sequences, source_db, 'tsv')
    
    # parse gain the gbff obtaining the features (tabular):
    feat_table = get_features(accession, gbff_file)

        
    # decide between 'locus_tag' and 'old_locus_tag' using intersection:
    i_locus     = set(feat_table['locus_tag'    ].to_list()).intersection(set(up_table['Gene Names (ordered locus)'].to_list()))
    i_old_locus = set(feat_table['old_locus_tag'].to_list()).intersection(set(up_table['Gene Names (ordered locus)'].to_list()))
    
    if len(i_old_locus) > len(i_locus):
        left_on = 'old_locus_tag'
    else: left_on = 'locus_tag'
    right_on = 'Gene Names (ordered locus)'
      
        
    if len(up_table)!=0 and (len(i_locus)>0 or len(i_old_locus)>0):  # at least one seq available in uniprot
    
        # merge based on different column names
        merge_table = pnd.merge(feat_table, up_table, how='inner', left_on=left_on, right_on=right_on)
        
        # select matching IDs:
        merge_table = merge_table[merge_table[left_on].notna()]
        # explanations on NCBI protein IDs: https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/reannotation/
        # if the .gbff is from GenBank, then 'protein_id' is the actual protein (example: CCC77580.1)
        # if the .gbff is from RefSeq, then 'protein_id' is a NR protein (prefix WP_)
        columns_to_keep = ['accession', 'seqid', 'gene_symbol', left_on, 'protein_id', 'Entry', 'RefSeq', 'KEGG', 'GeneID']
        merge_table = merge_table[columns_to_keep]
        
        
    else:  # uniprot data not available: 
        merge_table = feat_table.copy()
        
        columns_to_keep = ['accession', 'seqid', 'gene_symbol', left_on, 'protein_id', ]
        merge_table = merge_table[columns_to_keep]
        
        # add missing columns: 
        for i in ['Entry', 'RefSeq', 'KEGG', 'GeneID']:
            merge_table[i] = None
        
        
    
    # rename, reorder and select columns: 
    if source_db == 'genbank':
        renamed_cols = ['accession', 'seqid', 'gene_symbol', 'gid', 'ncbiprotein', 'uniprot', 'refseq', 'kegg', 'ncbigene']
    else:  # source_db == 'refseq'
        merge_table = merge_table.drop(columns='RefSeq')   # in this case it would be redundant
        # in case of 'refseq', the crossref 'ncbiprotein' is not available 
        merge_table['ncbiprotein'] = None
        renamed_cols = ['accession', 'seqid', 'gene_symbol', 'gid', 'refseq', 'uniprot', 'kegg', 'ncbigene', 'ncbiprotein']  
    # finally rename columns
    merge_table = merge_table.rename(columns=dict(zip(merge_table.columns, renamed_cols)))
    # reorder and select columns:
    merge_table = merge_table[['accession', 'gid', 'refseq', 'ncbiprotein', 'ncbigene', 'kegg', 'uniprot' ]]
    
    
    
    # if source_db=='genbank' , then the 'WP_' accessions of 'refseq' must be taken from 'up_table'.
    if source_db=='genbank':
        formatted_col = []
        for i in merge_table['refseq']:
            wp_code = None   # default
            if type(i) == str:  # avoid empty cells: 
                for code in i.split(';'):
                    #if code.startswith('WP_'):  # Memote is not yet ready for the new format.
                    if any([code.startswith(i) for i in "AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP".split('|')]):
                        wp_code = code
                        break
            formatted_col.append(wp_code)
        # replace with formatted version:
        merge_table['refseq'] = formatted_col
        
        
    # finally, the 'ncbigene' and 'kegg' columns need to be formatted:
    def simple_formatter(col):
        formatted_col = []
        for i in col: 
            formatted_cell = None
            if type(i) == str:  # avoid emtpy cells:
                formatted_cell = i.split(';')[0]
            formatted_col.append(formatted_cell)
        return formatted_col
    merge_table['ncbigene'] = simple_formatter(merge_table['ncbigene'])
    merge_table['kegg'] = simple_formatter(merge_table['kegg'])
    

    
    merge_table.to_csv(f'working/gannots/gannots_{accession}.csv')
    return merge_table



def create_proteome_from_gbff(accession, gbff_file, gannots, outfolder):
    
    
    # use SeqIO.read() if you expect only one record in the GenBank file.
    # alternatively, use SeqIO.parse() if there are multiple records in the file.
    records = SeqIO.parse(gbff_file, 'genbank')
                              

    sr_list = []
    for record in records:        
        
        # access features:
        for feature in record.features:
            if feature.type == 'CDS':
                try: location = feature.location
                except: location = None
                try: gene_symbol = feature.qualifiers['gene'][0]
                except: gene_symbol = None
                try: aa_seq = feature.qualifiers['translation'][0]
                except: 
                    try: aa_seq = str(feature.extract(record.seq).translate())  # reverse complement automatically handled
                    except: aa_seq = None
                try: locus_tag = feature.qualifiers['locus_tag'][0]
                except: locus_tag = None
                if locus_tag != None and (type(locus_tag)==int or str(locus_tag).isdigit()):  
                    # gene nume is just a number (the genbank may come from prodigal/bactabolize like pipelines)
                    # probably the gene ID will conflict with genes from other strains. 
                    # Therefore we prepend its accession:
                    locus_tag = accession + '_' + record.id + '_' + str(locus_tag)
                    locus_tag = locus_tag.replace('.', '').replace('-', '')  # remove undesired chars
                try: old_locus_tag = feature.qualifiers['old_locus_tag'][0]
                except: old_locus_tag = None
                try: protein_id = feature.qualifiers['protein_id'][0]
                except: protein_id = None
                
                
                # get the right 'gid'
                if old_locus_tag != None and old_locus_tag in gannots['gid'].to_list():
                    gid = old_locus_tag
                else: gid = locus_tag
                
                
                # append the seqrecord:
                if aa_seq != None and gid != None:
                    sr = SeqRecord.SeqRecord(Seq.Seq(aa_seq), id=gid, description=f'(file: {gbff_file} sequence: {record.id})')
                    sr_list.append(sr)
            
    
    # finally write the proteome in fasta format:
    proteome_filepath = outfolder + f'{accession}.faa'
    with open(proteome_filepath, 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
        
    return proteome_filepath
        
                
        
def create_proteome_and_gannots(gbff_file, outfolder):
    
    
    # get the "accession" for this file (just the simple filename)
    basename = os.path.basename(gbff_file)
    accession, _ = os.path.splitext(basename)
    
                       
    # create gene annotation pnd.DataFrame (gannots):
    gannots = get_gannots_tabular(accession, gbff_file, outfolder)
    
    
    # create proteome fasta file:
    proteome_filepath = create_proteome_from_gbff(accession, gbff_file, gannots, outfolder)
    
    
    return proteome_filepath, gannots



def check_cached_parsing(logger, species_to_genbank, outdir): 
    
    
    # verify the presence of needed files, and return 
    # the 'species_to_proteome' if they are all available.
    species_to_proteome = {}
    
    
    # if present , load the 'all_gannots' dataframe: 
    if not os.path.exists(outdir + 'gannots.csv'):
        return False
    # "low_memory=False" to prevent "DtypeWarning: Columns ... have mixed types."
    all_gannots = pnd.read_csv(outdir + 'gannots.csv', index_col=0, low_memory=False)
    avail_accessions = set(all_gannots['accession'].to_list())
        
    
    # check the presence of proteomes and gene annotation files: 
    required_accessions = set()
    for species in species_to_genbank.keys(): 
        species_to_proteome[species] = []
        
        for file in species_to_genbank[species]:
            basename = os.path.basename(file)
            accession, _ = os.path.splitext(basename)
            required_accessions.add(accession)
            
            # return empty if not found
            if not os.path.exists('working/proteomes/' + accession + '.faa'):
                return False
            if not os.path.exists('working/gannots/gannots_' + accession + '.csv'):
                return False
            
            # populate the 'species_to_proteome' dict if files are already available:
            species_to_proteome[species].append('working/proteomes/' + accession + '.faa')
            
            
    # the 'all_gannots' dataframe must have exaclty the same accessions, no more, no less:
    if required_accessions != avail_accessions:
        return False
            
            
    # save the dictionary to disk: 
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome, file)
           
            
    return True
            
    

def handle_manual_genbanks(logger, genbanks, outdir, metadata):
    
    
    # create a species-to-genome dictionary
    species_to_proteome = {}
    species_to_genbank = {}
    logger.debug(f"Checking the formatting of the provided -gb/--genbanks attribute...") 
    
    
    
    # PART A) populate the species_to_genbank dictionary.
    
    # check if the user specified a folder:
    if os.path.exists(genbanks):
        if os.path.isdir(genbanks):
            if genbanks[-1] != '/': genbanks = genbanks + '/'
            files = glob.glob(genbanks + '*')
            species_to_genbank['Spp'] = files
    
    elif '+' in genbanks and '@' in genbanks: 
        for species_block in genbanks.split('+'):
            species, files = species_block.split('@')
            for file in files.split(','): 
                if not os.path.exists(file):
                    logger.error("The following file provided in -gb/--genbanks does not exists: " + file)
                    return 1
            species_to_genbank[species] = files.split(',')
            
    else: # the user has just 1 species
        for file in genbanks.split(','): 
            if not os.path.exists(file):
                logger.error("The following file provided in -gb/--genbanks does not exists: " + file)
                return 1
        species_to_genbank['Spp'] = ganbanks.split(',')

    
    # report a summary of the parsing: 
    logger.info(f"Inputted {len(species_to_genbank.keys())} species with well-formatted paths to genbanks.") 
    
    

    # PART B) parse each genbank file, and populate the species_to_proteome dictionary.
    # Moreover, concat the gannots to create a single annotation table. 
    logger.info("Parsing genbanks file to extract proteomes and gene annotations...")
    os.makedirs('working/proteomes/', exist_ok=True)
    os.makedirs('working/gannots/', exist_ok=True)
    
    
    # check if files were already computed:
    if check_cached_parsing(logger, species_to_genbank, outdir):
        logger.info("Found all the needed files already computed. Skipping this step.")
        
        return 0
    
    
    # parse each genbank file:
    all_gannots = []
    for species in species_to_genbank.keys():
        created_files = []
        for file in species_to_genbank[species]:
            
            # all proteomes and gannots will be stored in 'working/proteomes/'.
            proteome_filepath, gannots = create_proteome_and_gannots(file, outfolder='working/proteomes/')
            created_files.append(proteome_filepath)
            all_gannots.append(gannots)
        species_to_proteome[species] = created_files
    logger.debug(f"Input genbanks converted to proteomes in ./working/proteomes/.")
    logger.debug(f"Created the species-to-proteome dictionary: {str(species_to_proteome)}.") 

    
    # create a single gene annotation table as final 'gempipe recon' output: 
    all_gannots = pnd.concat(all_gannots)
    all_gannots = all_gannots.astype(str)  # 'ncbigene' is still formatted as number
    all_gannots = all_gannots.replace('nan', None)
    #all_gannots['refseq'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in all_gannots['refseq'] ]
    all_gannots['ncbigene'] = [i.rsplit('.', 1)[0] if type(i)==str else None for i in all_gannots['ncbigene'] ]
        
    
    all_gannots.to_csv(outdir + 'gannots.csv')
    logger.debug(f"Saved the genes annotation table to {outdir + 'all_gannots.csv'}.")
    
    
    # save the dictionary to disk: 
    with open('working/proteomes/species_to_proteome.pickle', 'wb') as file:
        pickle.dump(species_to_proteome, file)
    logger.debug(f"Saved the species-to-proteome dictionary to file: ./working/proteomes/species_to_proteome.pickle.")
    
    
    # Create the genomes/genomes.csv like if genomes were downloaded from NCBI.
    # Useful during plot generation.
    # Warning: the same columns are used in get_metadata_table(). But here only 2 can be filled: 'organism_name' and 'strain_isolate'.
    get_genomes_csv(source='species_to_proteome')
    response = update_metadata_manual(logger, metadata, source='species_to_proteome')
    if response==1: return 1
    
    
    
    return 0