import json

import pandas as pnd
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap




def animatrix(
    tree_original='fastani/ANIclustermap_dendrogram.nwk', 
    triangular='fastani/ANIclustermap_matrix.tsv',
    type_strain=None, 
    legend_ratio=0.25, legend_title='species', genomes=None, cellannot=True, colorannot='niche', 
    replace0=80, fastmode=False, niche=False, excludeniche=[], nichewidth=2,
    outfile=None, verbose=False ):
    """Create a ANI plot starting from the ANIclustermap outputs.
    
    ANIclustermap can be found at github.com/moshi4/ANIclustermap.
    
    Args:
        tree_original (str): filepath to the `ANIclustermap_dendrogram.nwk` file produced by ANIclustermap.
        triangular (str): filepath to the `ANIclustermap_matrix.tsv` file produced by ANIclustermap.
        type_strain (str): accession of the genome from type strain to be used for a taxonoy-based filtering.
            All genomes having ANI < 95 respect to the type strain are discarded. 
        legend_ratio (float): space reserved for the legend.
        genomes (pandas.DataFrame): having at last the following columns `accession`, `strain`, `species`, `niche`. 
            The report produced by `gempipe recon` is fully compatible.
        cellannot (bool): if `True`, annotate each cell.
        colorannot (str): choose between 'species' and 'niche'.
        replace0 (float): replace 0s in the `ANIclustermap_matrix.tsv` with the provided value.
            This may improve the coloring. 
        fastmode (bool): if `True`, don't draw the dendrogram and make the picture smaller. 
        niche (bool): if `True`, draw squares according to the 'niche' attribute contained in 'genomes'.
        excludeniche (list): list of niches to exclude from the representation (when using 'niche').
            Bug: no more than 1 key is allowed.
        nichewidth (int): width of the border of squares (when using 'niche').
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages.

    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the ANI tree and associated heatmap.
    """
    

    # (1) load the tree produced by aniclustermap
    tree_original = Phylo.read(tree_original, "newick")
    # Sometimes, when used after 'ncbi_genome_download', leaves coud be formatted
    # like 'GCA_000010005.1_ASM1000v1_genomic', while the matrix like 'GCA_000010005.1'.
    # Here the formatting is corrected if needed:
    for leaf in tree_original.get_terminals():
        splitted = leaf.name.split('_')
        if len(splitted) > 2:  # (like in GCA_000010005.1)
            gca_gcf_index = None
            for i, split in enumerate(splitted): 
                if split in ['GCA', 'GCF']:
                    gca_gcf_index = i
            if gca_gcf_index != None:
                leaf.name = f'{splitted[gca_gcf_index]}_{splitted[gca_gcf_index +1]}'
    # get the leaves from top to bottom. It will be used to sort DataFrames later. 
    ord_leaves = [leaf.name for leaf in tree_original.get_terminals()]


    # (2) load the triangular matrix produced by aniclustermap
    triangular = pnd.read_csv(triangular, sep='\t')
    # improve columns and rows labeling: 
    columns = triangular.columns.to_list()
    columns = ['_'.join(c.split('_', 2)[:2]) for c in columns ]
    triangular.columns = columns
    triangular.index = columns
    
    
    # (2bis) filter based on taxonomy:
    # Given an type strain for a species. 
    if type_strain != None: 
        lost_acc = triangular[triangular[type_strain] < 95].index
        triangular = triangular[triangular[type_strain] >= 95]
        triangular = triangular[triangular.index]  # columns == index
        if verbose: 
            print(f"{len(lost_acc)} genomes had ANI < 95 respect to {type_strain}, thus were removed.")
        for leaf in lost_acc:
            tree_original.prune(leaf)
        ord_leaves = [leaf.name for leaf in tree_original.get_terminals()]
    if replace0 !=None: 
        triangular = triangular.replace(0, replace0)

 
    # (3) create the frame
    if genomes is not None:
        tree_ratio = 0.32
    else:
        tree_ratio = 0.19
        legend_ratio = 0
    niche_bar_ratio = 0.04 if niche else 0
    proportions = [tree_ratio, niche_bar_ratio, 0.02, 1.0,  0.02, 0.04, legend_ratio ]
    height = 0.2*len(ord_leaves)
    width = height * sum(proportions)
    if fastmode:
        # height : width = 10 : x 
        width = width * 10 / height
        height = 10
    fig, axs = plt.subplots(
        nrows=1, ncols=len(proportions),
        figsize=(width, height), # global dimensions.
        gridspec_kw={'width_ratios': proportions}) # suplots width proportions.
    # adjust the space between subplots: 
    plt.subplots_adjust(wspace=0, hspace=0)
    axs[1].axis('off')  # remove frame and axis
    axs[2].axis('off')  # remove frame and axis
    axs[4].axis('off')  # remove frame and axis
    
    
    # (4) get the colors and labels
    if genomes is not None:
        genomes = genomes.copy().set_index('accession', drop=False)
        ord_leaves_present = [i for i in ord_leaves if i in genomes.index]  # drop low-quality genomes
        if verbose: 
            print("Leaves accessions not found in 'genomes':", [i for i in ord_leaves if i not in genomes.index])
        genomes = genomes.loc[ord_leaves_present, ]  # drop low-quality genomes
        genomes['label'] = ''
        for accession, row in genomes.iterrows(): 
            genomes.loc[accession, 'label'] = f"{row['species']} {row['strain']} ({row['niche']})"
        if   colorannot=='species': color_key = 'species'
        elif colorannot=='niche': color_key = 'niche'
        else: 
            print("WARNING: wrong 'colorannot' parameter.")
            color_key = 'species'
        key_to_color = {key: f'C{number}' for number, key in enumerate(sorted(genomes.sort_index()[color_key].unique()))} 
        if  niche:
            for n in excludeniche: 
                key_to_color[n] = 'white'
        acc_to_color = genomes[color_key].map(key_to_color).to_dict()
    else:
        acc_to_color = None
            
        
    # (4) draw the dendrogram with colors and niche
    def get_label(leaf):
        if fastmode:
            return ''
        if leaf.name != None:
            if genomes is not None:
                row = genomes[genomes['accession']==leaf.name].iloc[0]
                return row['label']
            else:
                return leaf.name
        else: 
            return ''
    def get_color(leaf_name):
        if leaf_name != '':
            if acc_to_color is not None:
                row = genomes[genomes['label']==leaf_name].iloc[0]
                return acc_to_color[row['accession']]
            else:
                return 'black'
        else: 
            return 'black'
    Phylo.draw(
        tree=tree_original, 
        axes=axs[0],
        label_func=get_label,
        label_colors=get_color,
        do_show=False)
    axs[0].axis('off')  # remove frame and axis:
    # make the tree closer to the heatmap:
    if genomes is not None:
        spacer = len(max(genomes['label'].to_list())) * 14
    else:
        spacer = len(max(ord_leaves)) * 10
    if fastmode:
        spacer = 0
    axs[0].set_xbound(lower=0, upper= max(tree_original.depths().values()) + spacer )
    
    
    # (5) heatmap
    ordered_triangular = triangular.loc[ord_leaves]
    heatmap = axs[3].matshow(
        ordered_triangular,
        cmap=plt.cm.plasma, # colormap.
        #vmin=0, vmax=100, # define ranges for the colormap.
        aspect='auto') # 'auto': fit in the axes; 'equal': squared pixels
    axs[3].axis('off')  # remove frame and axis
    # draw the heatmap colormap (legend) in a separate ax
    plt.colorbar(heatmap, cax=axs[5]) 
    # annotate each cell if requested:
    if cellannot and not fastmode: 
        for i, row in enumerate(ordered_triangular.index):
            for j, col in enumerate(ordered_triangular.columns):
                value = ordered_triangular.loc[row, col]
                axs[3].text(j, i, str(round(value)), ha='center', va='center', color='black', fontsize=6)
    # annotate niches if requested: 
    if niche: 
        for n in genomes['niche'].unique():
            if n in excludeniche: 
                continue
            subset = genomes[genomes['niche']==n]
            first_acc = subset.iloc[0]['accession']
            last_acc = subset.iloc[-1]['accession']
            first_index = genomes.index.get_loc(first_acc) 
            last_index = genomes.index.get_loc(last_acc) 
            min_value = min([first_index, last_index])
            max_value = max([first_index, last_index])
            coords = (len(genomes)-min_value-(max_value-min_value)-1-0.5, min_value-0.5, max_value-min_value+1, max_value-min_value+1)
            rect = Rectangle((coords[0], coords[1]), coords[2], coords[3], linewidth=nichewidth, edgecolor=key_to_color[n], facecolor='none')
            axs[3].add_patch(rect)


    # (6) legend
    if genomes is not None and not niche and not fastmode:
        patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(sorted(genomes.sort_index()['species'].unique()))]
        l1 = plt.legend(handles=patches, title=legend_title, loc='center right')
        axs[6].add_artist(l1)  # l2 implicitly replaces l1
    if niche: 
        patches = [Patch(facecolor=key_to_color[n], label=n, ) for n in sorted(genomes['niche'].unique()) if n not in excludeniche]
        l2 = plt.legend(handles=patches, title=legend_title, loc='center right')
        axs[6].add_artist(l2)  # l2 implicitly replaces l1
    axs[6].axis('off')  # remove frame and axis
    
    
    # (7) create the colobar for 'niche'
    if niche: 
        colors_list = list(key_to_color.values())
        custom_cmap = LinearSegmentedColormap.from_list('CustomColormap', colors_list, N=256)
        matshow_group = [list(key_to_color.values()).index(acc_to_color[acc]) for acc in ord_leaves] 
        matshow_df = pnd.DataFrame({'accession': ord_leaves, 'group': matshow_group}).set_index('accession')
        clusters_matshow = axs[1].matshow(
            matshow_df[['group']],
            cmap= custom_cmap, 
            aspect='auto')
    
    
    # (7) save to disk
    # bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
    
    
    return fig, triangular



def phylogenomics(
    thirdparty_pam='roary/_1735552997/gene_presence_absence.csv', mode='roary', 
    newick='raxml_ng/core_gene_alignment.aln.raxml.bestTree', 
    legend_ratio=0.25, legend_title='species', genomes=None, showtiles=True, support_values=True, outgroup=None, 
    outfile=None, verbose=False ):
    """Create a phylogenomics tree starting from a pangenome analysis (e.g. the Roary outputs).
    
    Roary can be found at github.com/sanger-pathogens/Roary.
    
    Args:
        thirdparty_pam (str): filepath to the presence/absence matrix created by the pangenomics pipeline (e.g. `gene_presence_absence.csv` from Roary).
        mode (str): pangenomics pipeline used (Only "roary" is supported for the moment).
        legend_ratio (float): space reserved for the legend.
        legend_title (str): title for the species legend.
        genomes (pandas.DataFrame): having at last the following columns `accession`, `strain`, `species`, `niche`. 
            The report produced by `gempipe recon` is fully compatible.
        showtiles (bool): if `True`, include a graphical representation of the `thirdparty_pam`.
        support_values (bool): if `True`, indicate the support values.
        outgroup (str): improve representation of the outgroup (if present) in the tiles.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages.

    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the phylogenomics tree.
    """
    
    
    # (1) load the tree:
    newick = Phylo.read(newick, 'newick')
    # get the leaves from top to bottom. It will be used to sort DataFrames later. 
    ord_leaves = [leaf.name for leaf in newick.get_terminals()]
    
    
    # (2) load the pam
    if mode=='roary':
        pam = pnd.read_csv(thirdparty_pam, na_filter=False, low_memory=False)   # 'low_memory=False' as columns have mixed types.
        # binarize the matrix
        pam_binary = pam.iloc[:, 14:len(pam.columns)].apply(lambda x: x.map(lambda y: 0 if y=='' else 1))
    elif mode=='proteinortho': 
        pam = pnd.read_csv(thirdparty_pam, sep='\t', na_filter=False, low_memory=False)
        pam.columns = [i[:-len('.faa')] if i.endswith('.faa') else i for i in pam.columns  ]
        # binarize the matrix
        pam_binary = pam.iloc[:, 3:len(pam.columns)].apply(lambda x: x.map(lambda y: 0 if y=='*' else 1))
    elif mode=='orthofinder': 
        pam = pnd.read_csv(thirdparty_pam, sep='\t', na_filter=False, low_memory=False)
        # binarize the matrix
        pam_binary = pam.iloc[:, 1:len(pam.columns)].apply(lambda x: x.map(lambda y: 0 if y=='' else 1))

    
    #  order by sum of cols
    pam_binary = pam_binary.iloc[pam_binary.sum(axis=1).sort_values(ascending=False).index, :].reset_index(drop=True)
    # 'genomes' contains all the genomes, while the 'pam_binary' and 'newick' are made with quality-filtered genomes.
    if outgroup != None: # bring outgroup singleton to the tail: 
        # (df.drop('col_1', axis=1) == 0) creates a boolean DataFrame where each column (other than 'outgroup') is checked to be 0.
        # .all(axis=1) ensures that all the columns (other than 'outgroup') are 0 in each row.
        outgroup_singletons = pam_binary[(pam_binary[outgroup] == 1) & (pam_binary.drop(columns=[outgroup]) == 0).all(axis=1)]
        pam_binary = pam_binary.drop(index=outgroup_singletons.index)
        pam_binary = pnd.concat([pam_binary, outgroup_singletons])
    if verbose: 
        print('pam_binary.shape', pam_binary.shape)
        
        
    # (3) create the frame
    if genomes is not None:
        tree_ratio = 0.60
    else:
        tree_ratio = 0.50
        legend_ratio = 0
    if showtiles: 
        proportions = [tree_ratio, tree_ratio, legend_ratio ]
    else:
        proportions = [tree_ratio, legend_ratio ]
    height = 0.2*len(ord_leaves)
    fig, axs = plt.subplots(
        nrows=1, ncols=len(proportions),
        figsize=(height * sum(proportions), height), # global dimensions.
        gridspec_kw={'width_ratios': proportions}) # suplots width proportions.
    # adjust the space between subplots: 
    plt.subplots_adjust(wspace=0, hspace=0)
    
    
    # (4) get the colors
    if genomes is not None:
        genomes = genomes.copy().set_index('accession', drop=False)
        genomes = genomes.loc[ord_leaves, ]  # drop low-quality genomes
        genomes['label'] = ''
        for accession, row in genomes.iterrows(): 
            genomes.loc[accession, 'label'] = f"{row['strain']} ({row['niche']})"
        key_to_color = {key: f'C{number}' for number, key in enumerate(sorted(genomes.sort_index()['species'].unique()))}  
        acc_to_color = genomes['species'].map(key_to_color).to_dict()
    else:
        acc_to_color = None
        
        
    # (5) plot the tree
    # define labels and label colors:
    def get_leaf_label(leaf):
        if leaf.name != None:
            if genomes is not None:
                row = genomes[genomes['accession']==leaf.name].iloc[0]
                return row['label']
            else:
                return leaf.name
        else: 
            return ''
    def get_color(leaf_name):
        if leaf_name != '':
            if acc_to_color is not None:
                row = genomes[genomes['label']==leaf_name].iloc[0]
                return acc_to_color[row['accession']]
            else:
                return 'black'
        else: 
            return 'black'
    def get_node_label(node):
        if node.confidence != None:
            if node.confidence != 1:
                return node.confidence
            else:
                return None   # remove 1 to improve readibility
        else:
            return None
    Phylo.draw(
        tree=newick, 
        axes=axs[0],
        label_func=get_leaf_label,
        label_colors=get_color,
        branch_labels=get_node_label if support_values else None,
        show_confidence=support_values, 
        do_show=False)
    axs[0].axis('off')  # remove frame and axis:
    # make the tree closer to the heatmap:
    if genomes is not None:
        spacer = len(max(genomes['label'].to_list())) * 0.01
    else:
        spacer = len(max(ord_leaves)) * 0.013
    axs[0].set_xbound(lower=0, upper= max(newick.depths().values()) + spacer )
    
    
    # (6) plot the tiles
    if showtiles:
        matshow_df = pam_binary.T.loc[ord_leaves]
        axs[1].matshow(
            matshow_df +0.2, # transposed and reordered matrix. 
            cmap=plt.cm.Greys, 
            vmin=0, vmax=1, 
            aspect='auto', 
        )
        axs[1].axis('off')  # remove frame and axis
        
    
    # (7) legend
    ax_legend = axs[2] if showtiles else axs[1]
    if genomes is not None:
        patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(genomes.sort_index()['species'].unique())]
        l1 = plt.legend(handles=patches, title=legend_title, loc='center right')
        ax_legend.add_artist(l1)  # l2 implicitly replaces l1
    ax_legend.axis('off')  # remove frame and axis
    
    
    # (8) save to disk
    # bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
        
        
    return fig



def unwrap_ncbidataset(filepath='./assembly_data_report.jsonl'):
    """Extract metadata from the output file of "datasets download genome taxon" ('assembly_data_report.jsonl').
    
    NCBI dataset CLI can be found at: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/
    Citation: 10.1038/s41597-024-03571-y.
    
    Args:
        filepath (str): filepath to the 'assembly_data_report.jsonl' file.

    Returns:
        pandas.DataFrame: unwrapped metadata table. 
    """

    json_objects = []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                json_object = json.loads(line.strip()) 
                json_objects.append(json_object)
            except json.JSONDecodeError:
                print(f"Error decoding JSON on line: {line}")

    metadata = pnd.DataFrame.from_records(json_objects)

    # unwrap metadata: 
    metadata_unwrap = []
    for index, row in metadata.iterrows(): 
        row_dict = {}
        for key1 in metadata.columns:
            if type(row[key1]) != dict: 
                row_dict[f'{key1}'] = row[key1]
            else: 
                for key2 in row[key1].keys():
                    if type(row[key1][key2]) != dict: 
                        row_dict[f'{key1}_{key2}'] = row[key1][key2]
                    else: 
                        for key3 in row[key1][key2].keys():
                            if type(row[key1][key2][key3]) != dict: 
                                row_dict[f'{key1}_{key2}_{key3}'] = row[key1][key2][key3]
                            else:
                                for key4 in row[key1][key2][key3].keys():
                                    if type(row[key1][key2][key3][key4]) != dict: 
                                        row_dict[f'{key1}_{key2}_{key3}_{key4}'] = row[key1][key2][key3][key4]
                                    else:
                                        for key5 in row[key1][key2][key3][key4].keys():
                                            row_dict[f'{key1}_{key2}_{key3}_{key4}_{key5}'] = row[key1][key2][key3][key4][key5]
        metadata_unwrap.append(row_dict)

    metadata_unwrap = pnd.DataFrame.from_records(metadata_unwrap)
    metadata_unwrap  = metadata_unwrap.sort_values(by='accession')
    metadata_unwrap = metadata_unwrap.set_index('accession', drop=True)
    
    return metadata_unwrap



def get_filtering_summary(working_dir='gempipe/working/', thr_N50=50000, thr_nc=200, thr_bm=2, thr_bf=100, verbose=True):
    """
    Get a summary of the metrics used to filter genomes.
    
    Args:
        working_dir (str): path to the gempipe's working directory.
        thr_N50 (int): N50 threshold.
        thr_nc (int): number of contigs threshold.
        thr_bm (int): BUSCO M% (missing) threshold.
        thr_bf (int): BUSCO F% (fragmented) threshold.

    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: metrics for all input genomes.
            - pandas.DataFrame: metrics only for retained genomes.
            - matplotlib.figure.Figure: histograms of metrics.
    """
    
    # load tables by reading the working directory: 
    genomes = pnd.read_csv(f'{working_dir}/genomes/genomes.csv', index_col=0)
    genomes = genomes.sort_values(by='assembly_accession')
    t_metrics = pnd.read_csv(f'{working_dir}/filtering/tmetrics.csv', index_col=0)
    t_metrics = t_metrics.set_index('accession', drop=True)
    b_metrics = pnd.read_csv(f'{working_dir}/filtering/bmetrics.csv', index_col=0)
    b_metrics = b_metrics.set_index('accession', drop=True)
    
    
    # create a table whowing all metrics ('summary_table')
    summary_table = pnd.concat([genomes.set_index('assembly_accession', drop=True), t_metrics], axis=1)
    summary_table = pnd.concat([summary_table, b_metrics], axis=1)
    summary_table['sum_len'] = round(summary_table['sum_len'] /  1000 / 1000, 3)
    summary_table = summary_table[['strain_isolate', 'organism_name', 'ncontigs', 'sum_len', 'N50', 'GC(%)', 'F', 'M']]
    summary_table = summary_table.rename(columns={'F': 'BUSCO_F%', 'M': 'BUSCO_M%'})
    

    # apply thresholds:
    summary_table_filt = summary_table[
        (summary_table['N50'] >= thr_N50) & \
        (summary_table['ncontigs'] <= thr_nc) & \
        (summary_table['BUSCO_M%'] <= thr_bm) & \
        (summary_table['BUSCO_F%'] <= thr_bf)]
    
    
    if verbose: 
        print(f"Removed after 'N50'>={thr_N50}: {len(summary_table[(summary_table['N50']<thr_N50)])} / {len(summary_table)} strains")
        print(f"Removed after 'ncontigs%'<={thr_nc}: {len(summary_table[(summary_table['ncontigs']>thr_nc)])} / {len(summary_table)} strains")
        print(f"Removed after 'BUSCO_M%'<={thr_bm}: {len(summary_table[(summary_table['BUSCO_M%']>thr_bm)])} / {len(summary_table)} strains")
        print(f"Removed after 'BUSCO_F%'<={thr_bf}: {len(summary_table[(summary_table['BUSCO_F%']>thr_bf)])} / {len(summary_table)} strains")
        
        print(f"Remaining: {len(summary_table_filt)} / {len(summary_table)} strains")
    
    
    # define the canvas:
    nrows = 4
    props = [1 for i in range(nrows)]
    fig, axs = plt.subplots(
        nrows=4, ncols=1,
        figsize=(8,6), # global dimensions.
        gridspec_kw={'height_ratios': props}) # suplots width proportions.
    # adjust the space between subplots: 
    plt.subplots_adjust(wspace=0, hspace=1)
    fig.supylabel('absolute counts')
    
    for i, metric, thr in zip(list(range(nrows)), ['N50','ncontigs', 'BUSCO_M%','BUSCO_F%'], [thr_N50, thr_nc, thr_bm, thr_bf]): 
        
        summary_table[metric].plot.hist(bins=200, ax=axs[i], ylabel='')
        if not (metric in ['BUSCO_M%','BUSCO_F%'] and thr==100): 
            axs[i].axvline(x=thr, color='red', linestyle='-', linewidth=0.5)
        axs[i].set_xlabel(metric)
        axs[i].set_facecolor('ghostwhite')
        if metric == 'N50':   # the only case where higher values are better
            axs[i].invert_xaxis()
    
    
    return summary_table, summary_table_filt, fig