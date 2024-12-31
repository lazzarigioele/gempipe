import pandas as pnd
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def aniclustermap_enhanced(
    tree_original='fastani/ANIclustermap_dendrogram.nwk', 
    triangular='fastani/ANIclustermap_matrix.tsv', 
    legend_ratio=0.25, genomes=None, cellannot=True, 
    outfile=None, verbose=False 
):
    """Create a ANI plot starting from the ANIclustermap outputs.
    
    ANIclustermap can be found at github.com/moshi4/ANIclustermap.
    
    Args:
        tree_original (str): filepath to the `ANIclustermap_dendrogram.nwk` file produced by ANIclustermap.
        triangular (str): filepath to the `ANIclustermap_matrix.tsv` file produced by ANIclustermap.
        legend_ratio (float): space reserved for the legend.
        genomes (pandas.DataFrame): having at last the following columns `assembly_accession`, `strain_isolate`, `organism_name`, `niche`. 
            The teble produced in `working/genomes/genomes.csv` is fully compatible.
        cellannot (bool): if `True`, annotate each cell.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages.

    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the ANI tree and associated heatmap.
    """
    

    # (1) load the tree produced by aniclustermap
    tree_original = Phylo.read(tree_original, "newick")
    # get the leaves from top to bottom. It will be used to sort DataFrames later. 
    ord_leaves = [leaf.name for leaf in tree_original.get_terminals()]


    # (2) load the triangular matrix produced by aniclustermap
    triangular = pnd.read_csv(triangular, sep='\t')
    # improve columns and rows labeling: 
    columns = triangular.columns.to_list()
    columns = ['_'.join(c.split('_', 2)[:2]) for c in columns ]
    triangular.columns = columns
    triangular.index = columns

    
    # (3) create the frame
    if genomes is not None:
        tree_ratio = 0.32
    else:
        tree_ratio = 0.19
        legend_ratio = 0
    proportions = [tree_ratio, 1.0,  0.02, 0.04, legend_ratio ]
    height = 0.2*len(ord_leaves)
    fig, axs = plt.subplots(
        nrows=1, ncols=5,
        figsize=(height * sum(proportions), height), # global dimensions.
        gridspec_kw={'width_ratios': proportions}) # suplots width proportions.
    # adjust the space between subplots: 
    plt.subplots_adjust(wspace=0, hspace=0)
    axs[2].axis('off')  # remove frame and axis
    
    
    # (4) get the colors
    if genomes is not None:
        genomes = genomes.copy().set_index('assembly_accession', drop=False)
        genomes = genomes.loc[ord_leaves, ]  # drop low-qiality genomes
        genomes['label'] = ''
        for accession, row in genomes.iterrows(): 
            genomes.loc[accession, 'label'] = f"{row['strain_isolate']} ({row['niche']})"
        key2color = {key: f'C{number}' for number, key in enumerate(genomes.sort_index()['organism_name'].unique())}   
        acc_to_color = genomes['organism_name'].map(key2color).to_dict()
    else:
        acc_to_color = None

        
    # (4) draw the dendrogram with colors and niche
    def get_label(leaf):
        if leaf.name != None:
            if genomes is not None:
                row = genomes[genomes['assembly_accession']==leaf.name].iloc[0]
                return row['label']
            else:
                return leaf.name
        else: 
            return ''
    def get_color(leaf_name):
        if leaf_name != '':
            if acc_to_color is not None:
                row = genomes[genomes['label']==leaf_name].iloc[0]
                return acc_to_color[row['assembly_accession']]
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
    axs[0].set_xbound(lower=0, upper= max(tree_original.depths().values()) + spacer )
    
    
    # (5) heatmap
    ordered_triangular = triangular.loc[ord_leaves]
    heatmap = axs[1].matshow(
        ordered_triangular,
        cmap=plt.cm.plasma, # colormap.
        #vmin=0, vmax=100, # define ranges for the colormap.
        aspect='auto') # 'auto': fit in the axes; 'equal': squared pixels
    axs[1].axis('off')  # remove frame and axis
    # draw the heatmap colormap (legend) in a separate ax
    plt.colorbar(heatmap, cax=axs[3]) 
    if cellannot: # annotate each cell:
        for i, row in enumerate(ordered_triangular.index):
            for j, col in enumerate(ordered_triangular.columns):
                value = ordered_triangular.loc[row, col]
                axs[1].text(j, i, str(round(value)), ha='center', va='center', color='black', fontsize=6)


    # (6) legend
    if genomes is not None:
        patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(genomes.sort_index()['organism_name'].unique())]
        l1 = plt.legend(handles=patches, title='submitted species', loc='center right')
        axs[4].add_artist(l1)  # l2 implicitly replaces l1
    axs[4].axis('off')  # remove frame and axis
    
    
    # (7) save to disk
    # bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
    
    
    return fig
