import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram, leaves_list
from sklearn.metrics import silhouette_score, silhouette_samples
import numpy as np    



def merge_tables(dict_tables):
    
    
    if isinstance(dict_tables, pnd.DataFrame):
        dict_tables = {'unknown_layer': dict_tables}
        
    
    # concat tables:
    tables_list = []
    for key, value in dict_tables.items(): 
        tables_list.append(value)
        
    # verify shape:
    for table in tables_list: 
        if set(list(tables_list[0].columns)) != set(list(table.columns)):
            print("ERROR: provided tables have different accessions.")
            return
      
    # concat the tables:
    data = pnd.concat(tables_list)
    
    # remove any ro conntaining NA in at least 1 column :
    rows_with_missing = list(data[data.isna().any(axis=1)].index)
    if rows_with_missing != []:
        print(f"WARNING: removing rows with missing values: {rows_with_missing}.")
    data = data.dropna()
    
    # verify binary format:
    binary = data.isin([0, 1]).all().all() or data.isin([False, True]).all().all()
    if not binary:
        print("ERROR: provided data are not binary (must be all {0,1} or all {False,True}.")
        return
    
    # start with all cells 0/1
    data.astype(int)
    
    # accessions as rows, features as columns: 
    data = data.T
    
    
    # covert to multi-layer dataframe: 
    for col in data.columns:
        for i, (name, table) in enumerate(dict_tables.items()):
            if col in table.index:
                data[col] = data[col].replace(1, i+1)  
                
                
    return data, dict_tables



def make_dendrogram(ax, dendrogram_data):
        
    # plot the dendrogram
    dn = dendrogram(
        dendrogram_data, ax=ax,
        orientation='left',
        color_threshold=0,
        above_threshold_color='black')


    # remove frame borders: 
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # remove ticks and markers: 
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)



def make_colorbar_clusters(ax, index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix):

    if acc_to_cluster!=None and cluster_to_color!=None: 
        
        # create the colors: 
        colors_list = list(cluster_to_color.values())
        custom_cmap = LinearSegmentedColormap.from_list('CustomColormap', colors_list, N=256)

        # create the dataframe: 
        ord_leaves = leaves_list(linkage_matrix)
        matshow_acc   = [index_to_acc[index] for index in ord_leaves]
        matshow_group = [acc_to_cluster[index_to_acc[index]] for index in ord_leaves]
        mathow_df = pnd.DataFrame({'accession': matshow_acc, 'group': matshow_group}).set_index('accession')
        mathow_df = mathow_df[::-1]  # leaves are drawn from bottom to top

        clusters_matshow = ax.matshow(
            mathow_df[['group']],
            cmap= custom_cmap, 
            aspect='auto')

    ax.axis('off')  # remove frame and axis
    
    
    
def make_colorbar_metadata(ax, derive_report, report_key, index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix):

    
    if isinstance(derive_report, pnd.DataFrame):
        
        if report_key not in derive_report.columns:
            print("WARNING: provided 'report_key' not found in 'derive_report' columns.")
            ax.axis('off')  # remove frame and axis
            return
        
        # define accession-to-colors:
        acc_to_colors = derive_report[report_key].map({species: f'C{number}' for number, species in enumerate(derive_report[report_key].unique())}).to_dict()    
        # create the colors: 
        colors_list = list(acc_to_colors.values())
        custom_cmap = LinearSegmentedColormap.from_list('CustomColormap', colors_list, N=256)

        # create the dataframe: 
        ord_leaves = leaves_list(linkage_matrix)
        matshow_acc   = [index_to_acc[index] for index in ord_leaves]
        matshow_group = [acc_to_cluster[index_to_acc[index]] for index in ord_leaves]
        mathow_df = pnd.DataFrame({'accession': matshow_acc, 'group': matshow_group}).set_index('accession')
        mathow_df = mathow_df[::-1]  # leaves are drawn from bottom to top

        clusters_matshow = ax.matshow(
            mathow_df[['group']],
            cmap= custom_cmap, 
            aspect='auto')

    ax.axis('off')  # remove frame and axis
    


def make_legends(ax, derive_report, report_key, cluster_to_color, dict_tables):
    
    # l1: species / niche
    if isinstance(derive_report, pnd.DataFrame):
        patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(derive_report[report_key].unique())]
        l1 = plt.legend(handles=patches, title=report_key, loc='upper left')  # , bbox_to_anchor=(1.05, 0.5)
        ax.add_artist(l1)  # l2 implicitly replaces l1
        
    
    # l2: clusters
    patches = [Patch(facecolor=color, label=f"Cluster_{cluster}", ) for cluster, color in cluster_to_color.items()]
    l2 = plt.legend(handles=patches, title='clusters', loc='center left')  # , bbox_to_anchor=(1.05, 0.5)
    ax.add_artist(l2)  # l2 implicitly replaces l1
    
    
    # l3: features
    if dict_tables != None:
        n_colors = len(list(dict_tables.keys()))  +1  # +1 for 'absence'
        viridis_discrete = plt.cm.get_cmap('viridis', n_colors) 
        viridis_discrete_rgb = viridis_discrete([i for i in range(n_colors)])
        patches = [Patch(facecolor=viridis_discrete_rgb[i+1], label=key) for i, key in enumerate(dict_tables.keys())]
        patches = [Patch(facecolor=viridis_discrete_rgb[0], label='absence')] + patches
        l3 = plt.legend(handles=patches, title='features', loc='lower left')  # , bbox_to_anchor=(1.05, 0.5)
        ax.add_artist(l3)  # l2 implicitly replaces l1
    
    ax.axis('off')  # remove frame and axis
    
    