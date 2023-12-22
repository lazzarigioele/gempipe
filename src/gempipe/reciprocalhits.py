

def task_brh(genome, args):
    
    # retrive the arguments:
    #pam = args['pam']
    #sequences_df = args['sequences_df']
    #acc_to_suffix = args['acc_to_suffix']
    #cluster_to_rep = args['cluster_to_rep']
    
    
    # get the basename without extension:
    basename = os.path.basename(genome)
    accession, _ = os.path.splitext(basename)
    
    
    # create subdir without overwriting: 
    os.makedirs(f'manual_brh/{accession}/', exist_ok=True)
            
    
    # create blast databases: 

    os.makedirs(f'manual_brh/{accession}/dbs/reference/', exist_ok=True)
    # shutil.copy is specified to copy permission bits, while shutil.copyfile is just to copy the file contents.
    shutil.copyfile(f'from_mendoza/protein_fasta.faa', f'manual_brh/{accession}/dbs/reference/protein_fasta.faa')  
    # '-parse_seqids' is required for later use of 'blastdbcmd'.
    command = f"""makeblastdb -in manual_brh/{accession}/dbs/reference/protein_fasta.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()

    os.makedirs(f'manual_brh/{accession}/dbs/{accession}/', exist_ok=True)
    # shutil.copy is specified to copy permission bits, while shutil.copyfile is just to copy the file contents.
    shutil.copyfile(f'prodigal/prokka/{accession}/{accession}.faa', f'manual_brh/{accession}/dbs/{accession}/{accession}.faa')  
    # '-parse_seqids' is required for later use of 'blastdbcmd'.
    command = f"""makeblastdb -in manual_brh/{accession}/dbs/{accession}/{accession}.faa -dbtype prot"""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()

    # STEP 2: perform 2 blastp

    command = f'''blastp \
        -query manual_brh/{accession}/dbs/reference/protein_fasta.faa \
        -db manual_brh/{accession}/dbs/{accession}/{accession}.faa \
        -out manual_brh/{accession}/align_ref_vs_acc.tsv \
        -outfmt "6 qseqid sseqid pident ppos length sstart send qcovs evalue bitscore"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()

    command = f'''blastp \
        -query manual_brh/{accession}/dbs/{accession}/{accession}.faa \
        -db manual_brh/{accession}/dbs/reference/protein_fasta.faa \
        -out manual_brh/{accession}/align_acc_vs_ref.tsv \
        -outfmt "6 qseqid sseqid pident ppos length sstart send qcovs evalue bitscore"
    '''
    process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()

    shutil.rmtree(f'manual_brh/{accession}/dbs/') # erease to save space. 


    header = "qseqid sseqid pident ppos length sstart send qcovs evalue bitscore".split(' ')
    align_ref_vs_acc = pnd.read_csv(f'manual_brh/{accession}/align_ref_vs_acc.tsv', names=header, sep='\t')
    align_acc_vs_ref = pnd.read_csv(f'manual_brh/{accession}/align_acc_vs_ref.tsv', names=header, sep='\t')


    df = []
    for cds in align_acc_vs_ref['qseqid'].unique():

        curr_hsps = align_acc_vs_ref[align_acc_vs_ref['qseqid'] == cds]
        curr_hsps_filt = curr_hsps[(curr_hsps['qcovs'] >= 70) & (curr_hsps['evalue'] <= 1e-5)]
        curr_hsps_filt_sort = curr_hsps_filt.sort_values(by='evalue', ascending=True)
        curr_hsps_filt_sort.reset_index(inplace=True, drop=True)
        try: best_hit = curr_hsps_filt_sort.loc[0, 'sseqid']
        except: 
            df.append({'cds': cds, 'ref': None, 'reciprocal': 'NA'})
            continue


        curr_hsps2 = align_ref_vs_acc[align_ref_vs_acc['qseqid'] == best_hit]
        curr_hsps_filt2 = curr_hsps2[(curr_hsps2['qcovs'] >= 70) & (curr_hsps2['evalue'] <= 1e-5)]
        curr_hsps_filt_sort2 = curr_hsps_filt2.sort_values(by='evalue', ascending=True)
        curr_hsps_filt_sort2.reset_index(inplace=True, drop=True)
        try: best_hit2 = curr_hsps_filt_sort2.loc[0, 'sseqid']
        except: 
            df.append({'cds': cds, 'ref': best_hit, 'reciprocal': 'NA'})
            continue


        if cds == best_hit2: df.append({'cds': cds, 'ref': best_hit, 'reciprocal': '<=>'})
        else: df.append({'cds': cds, 'ref': best_hit, 'reciprocal': '=>'})



    df = pnd.DataFrame.from_records(df)
    df.to_csv(f'manual_brh/{accession}/results.csv', sep='\t')
    



def perform_brh(logger, cores, ref_proteome): 
    
    
    # some log messages:
    logger.info("Performing the best reciprocal hits (BRH) alignment against the reference proteome...")
    
    
    # create sub-directories without overwriting:
    os.makedirs('working/brh/', exist_ok=True)

    
    # TODO
    """
    # check if it's everything pre-computed
    response = check_cached(
        logger, pam_path='working/rec_overlap/pam.csv',
        summary_path='working/rec_overlap/summary.csv',
        imp_files = [
            'working/rec_overlap/sequences.csv',
            'working/rec_overlap/seq_to_coords.pickle',])
    if response == 0: 
        return 0
    """
    
    
    # TODO
    """
    # load the assets to form the args dictionary:
    pam = pnd.read_csv('working/rec_masking/pam.csv', index_col=0)
    sequences_df = pnd.read_csv('working/rec_masking/sequences.csv', index_col=0)
    with open('working/clustering/acc_to_suffix.pickle', 'rb') as handler:
        acc_to_suffix = pickle.load(handler)
    with open('working/clustering/cluster_to_rep.pickle', 'rb') as handler:
        cluster_to_rep = pickle.load(handler)
    """
        
        
    # load the previously created species_to_genome: 
    with open('working/genomes/species_to_genome.pickle', 'rb') as handler:
        species_to_genome = pickle.load(handler)
    
    
    # create items for parallelization: 
    items = []
    for species in species_to_genome.keys(): 
        for genome in species_to_genome[species]: 
            items.append(genome)
    

    # randomize and divide in chunks: 
    chunks = chunkize_items(items, cores)
    
    
    # initialize the globalpool:
    globalpool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    
    
    # start the multiprocessing: 
    results = globalpool.imap(
        load_the_worker, 
        zip(chunks, 
            range(cores), 
            itertools.repeat(['cds', 'accession', 'aaseq']), 
            itertools.repeat('cds'), 
            itertools.repeat(logger), 
            itertools.repeat(task_brh),  # ...TODO
            itertools.repeat({'pam': pam, 'sequences_df': sequences_df, 'acc_to_suffix': acc_to_suffix, 'cluster_to_rep': cluster_to_rep}),
        ), chunksize = 1)
    all_df_combined = gather_results(results)
    
    
    # TODO
    # save tabular results:
    all_df_combined
    
    
    # empty the globalpool
    globalpool.close() # prevent the addition of new tasks.
    globalpool.join() 