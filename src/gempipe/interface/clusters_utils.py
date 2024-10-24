import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram, leaves_list
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
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



def sort_by_leaves(data, linkage_matrix, index_to_acc):
    
    # How to get the leaves order: 
    ord_leaves = leaves_list(linkage_matrix)
    ord_leaves = np.flip(ord_leaves)  # because leaves are returned in the inverse sense.
    ord_leaves = [index_to_acc[i] for i in ord_leaves]  # convert index as number to index as accession
    ord_data = data.loc[ord_leaves, :]  # reordered dataframe.
    
    return ord_data



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



def make_colorbar_clusters(ax, ord_data, acc_to_cluster, cluster_to_color):

    if acc_to_cluster!=None and cluster_to_color!=None: 
        
        # create the colors: 
        colors_list = list(cluster_to_color.values())
        custom_cmap = LinearSegmentedColormap.from_list('CustomColormap', colors_list, N=256)

        
        # create a dataframe (matshow_df) with a single column ('group'):
        matshow_acc = [acc for acc in ord_data.index]
        matshow_group = [acc_to_cluster[acc] for acc in ord_data.index]   
        matshow_df = pnd.DataFrame({'accession': matshow_acc, 'group': matshow_group}).set_index('accession')
        
        clusters_matshow = ax.matshow(
            matshow_df[['group']],
            cmap= custom_cmap, 
            aspect='auto')

    ax.axis('off')  # remove frame and axis
    
    
    
def make_colorbar_metadata(ax, ord_data, derive_report, report_key):

    
    if isinstance(derive_report, pnd.DataFrame):
        
        if report_key not in derive_report.columns:
            print("WARNING: provided 'report_key' not found in 'derive_report' columns.")
            ax.axis('off')  # remove frame and axis
            return
        
        # define accession-to-colors:
        key_to_color = {key: f'C{number}' for number, key in enumerate(derive_report[report_key].unique())}   # 'key' is eg 'species'.
        acc_to_color = derive_report[report_key].map(key_to_color).to_dict()    
        # create the colors: 
        colors_list = list(key_to_color.values())
        custom_cmap = LinearSegmentedColormap.from_list('CustomColormap', colors_list, N=256)

      
        # create a dataframe (matshow_df) with a single column ('group'):
        matshow_acc = [acc for acc in ord_data.index]
        matshow_group = [int(acc_to_color[acc][1:]) for acc in ord_data.index]   # int(cell[1:]) convert color for example rom 'C1' to 1.
        matshow_df = pnd.DataFrame({'accession': matshow_acc, 'group': matshow_group}).set_index('accession')
                
        clusters_matshow = ax.matshow(
            matshow_df[['group']],
            cmap= custom_cmap, 
            aspect='auto')

    ax.axis('off')  # remove frame and axis
    


def make_legends(ax, derive_report, report_key, cluster_to_color, dict_tables):
    
    # l1: species / niche
    if isinstance(derive_report, pnd.DataFrame):
        patches = [Patch(facecolor=f'C{number}', label=species, ) for number, species in enumerate(derive_report[report_key].unique())]
        l1 = plt.legend(handles=patches, title=report_key, loc='upper left')  
        ax.add_artist(l1)  # l2 implicitly replaces l1
        
    
    # l2: clusters
    if cluster_to_color != None: 
        patches = [Patch(facecolor=color, label=f"Cluster_{cluster}", ) for cluster, color in cluster_to_color.items()]
        l2 = plt.legend(handles=patches, title='clusters', loc='center left')
        ax.add_artist(l2)  # l2 implicitly replaces l1
    
    
    # l3: features
    if dict_tables != None:
        n_colors = len(list(dict_tables.keys()))  +1  # +1 for 'absence'
        viridis_discrete = plt.cm.get_cmap('viridis', n_colors) 
        viridis_discrete_rgb = viridis_discrete([i for i in range(n_colors)])
        patches = [Patch(facecolor=viridis_discrete_rgb[i+1], label=key) for i, key in enumerate(dict_tables.keys())]
        patches = [Patch(facecolor=viridis_discrete_rgb[0], label='absence')] + patches
        l3 = plt.legend(handles=patches, title='features', loc='lower left')  
        ax.add_artist(l3)  # l2 implicitly replaces l1
    
    ax.axis('off')  # remove frame and axis
    
    
    
def subset_k_best(data_bool, k, acc_to_cluster, derive_report, report_key):
    
    if acc_to_cluster==None and not isinstance(derive_report, pnd.DataFrame):
        print("WARNING: to use 'k', 'acc_to_cluster' or 'derive_report' must be specified.")
        return data_bool
    
    
    # constant columns have already been dropped
    # The classical 'X' (features) and 'y' (classification) datasets are prepared. 
    X = data_bool
    if   isinstance(derive_report, pnd.DataFrame):  # priority
        y = [derive_report.loc[acc, report_key] for acc in data_bool.index]
        print(f"WARNING: feature selection to distinguish between '{report_key}'.")
    elif acc_to_cluster != None:
        y = [f'Cluster_{acc_to_cluster[acc]}' for acc in data_bool.index]

    # with SelectKBet, we can use different scoring functiond. 
    # The 'f_classif' (which corresponds to an ANOVA F-value or ANOVA F-statistic) works well for binary classification.
    # However, 'f_classif' is not ideal for multiclass classification.
    # The 'mutual_info_classif' should work better for multiclass problems. 
    selector = SelectKBest(score_func=mutual_info_classif, k=k)

    # Fit and transform the data to keep only the K best features.
    # Possible warnings could be:
    # /opt/conda/lib/python3.8/site-packages/sklearn/feature_selection/_univariate_selection.py:113: 
    # RuntimeWarning: invalid value encountered in divide f = msb / msw
    # This warning is raised because some of the data contains zero variance or constant values, 
    # leading to a division by zero or other numerical issues. 
    X_new = selector.fit_transform(X, y)  # still not used

    
    # The selected features can be accessed using selector.get_support()
    sel_feat_indices = selector.get_support(indices=True)  # <class 'numpy.ndarray'>
    X_sub = X.iloc[:, list(sel_feat_indices)]
    
    
    return X_sub
    
        