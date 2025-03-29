import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram, leaves_list
from sklearn.metrics import silhouette_score, silhouette_samples
import numpy as np    


from gempipe.interface.clusters_utils import merge_tables, sort_by_leaves, make_dendrogram, make_colorbar_clusters, make_colorbar_metadata, make_legends, subset_k_best

    
    
def silhouette_analysis(
    tables, figsize = (10,5), drop_const=True, ctotest=None, forcen=None, 
    derive_report=None, report_key='species', excludekeys=[],
    legend_ratio=0.7, outfile=None, verbose=False, anchor=[None, None, None], key_to_color=None):
    """Perform a silhuette analysis to detect the optimal number of clusters. 
    
    Args:
        tables (pnd.DataFrame): feature tables with genome accessions are in columns and features are in rows. 
            Can also be a dictionary of feature tables (example: ``{'auxotrophies': aux_df, 'substrates': sub_df})``. 
            In this case, any number of tables (pandas.DataFrame) can be used. 
            For each table, genome accessions are in columns, features are in rows.
            Directly compatible tables are: `rpam.csv`, `cnps.csv`, and `aux.csv` (all produced by `gempipe derive`).
        figsize (int, int): width and height of the figure.
        drop_const (bool): if `True`, remove constant features.
        ctotest (list): number of clusters to test (example: ``[5,7,10]`` to test five, seven and ten clusters).
            If `None`, all the combinations from 2 to the number of accessions -1 will be used.
        forcen (int): force the number of cluster, otherwise the optimal number will picked up according to the sihouette value. 
        derive_report (pandas.DataFrame): report table for the generation of strain-specific GSMMs, made by `gempipe derive` in the output directory (`derive_strains.csv`). 
        excludekeys (list): keys (iches/species) not to show in legend. 
            Bug: no more than 1 key is allowed.
        report_key (str): name of the attribute (column) appearing in `derive_report`, to be compared to the metabolilc clusters.
            Usually it is 'species' or 'niche'.
        legend_ratio (float): space reserved for the legend.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages.
        anchor (list): list of tuples (X,Y) for customixing the position of legends. 
            ``None`` will leave default positioning.
        key_to_color (dict): dict mapping each category in `report_key` to a color in the format ([0:1],[0:1],[0:1]).
            ``None`` will leave default color and order in the legend. 
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the sinhouette analysis.
            - dict: genome-to-cluster associations.
            - dict: an RGB color for each cluster.
    """
    
    
    def create_silhuette_frame(figsize):
    
        # create the subplots: 
        fig, axs = plt.subplots(
            nrows=1, ncols=10, 
            figsize=figsize, # global dimensions.  
            gridspec_kw={'width_ratios': [0.46, 0.02, 0.46, 0.02, 0.3, 0.04, 0.02, 0.04, 0.02, legend_ratio]}) # suplots width proportions. 
        # adjust the space between subplots: 
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[1].axis('off')  # remove frame and axis
        axs[3].axis('off')  # remove frame and axis
        axs[6].axis('off')  # remove frame and axis
        axs[8].axis('off')  # remove frame and axis

        return fig, axs



    def make_plot_ncluster_comparison(ax, num_clusters_vector, silhouette_avg_scores, opt_n_clusters, forcen, verbose):

        # Plot the silhouette scores against the number of clusters (threshold values)
        ax.plot(num_clusters_vector, silhouette_avg_scores, marker='o')
        ax.set_xlabel('N clusters')
        ax.set_ylabel('Average Silhouette Score')
        ax.grid(True)

        if verbose: print(f"Optimal number of clusters: {opt_n_clusters}")
        ax.axvline(x=opt_n_clusters if forcen==None else forcen, color='red', linestyle='--')



    def make_plot_silhouette_coeff(ax, opt_n_clusters, silhouette_scores, clusters):

        # Given a fixed number of cluster (ie the optimal number of clusters),
        # extract the datapoint belonging to each of the clusters and show its associates silhouette score. 
        y_lower = 0
        cluster_to_color = {}
        for i in range(opt_n_clusters):

            # scores of the datapoint inside the cluster i: 
            cluster_i_scores = silhouette_scores[clusters == i]
            cluster_i_scores.sort()  # sort from smallest to biggest 

            size_cluster_i = len(cluster_i_scores)
            # get the limits for this polygon:
            y_upper = y_lower + size_cluster_i
            # get the color for this cluster/polygon: 
            color = plt.cm.Spectral(float(i) / opt_n_clusters)
            cluster_to_color[i+1] = color

            ax.fill_betweenx(np.arange(y_lower, y_upper),
                              0, cluster_i_scores,
                              facecolor=color, edgecolor=color, alpha=1.0)
            ax.text(0, (y_lower + y_upper -1)/2, f'Cluster_{i+1}', va='center', ha='left')
            y_lower = y_upper + -1  # no space between clusters

        ax.set_xlabel('Silhouette Coefficient')
        ax.set_ylabel('')
        ax.set_yticks([]) 
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        #ax.set_title('Silhouette Plot for {} Clusters'.format(opt_n_clusters))
        ax.set_facecolor('#f8f8f8')  # Light gray background

        return cluster_to_color



    # START:
    # format input tables:
    data, dict_tables = merge_tables(tables)
    data_bool = data.astype(bool)   # convert multi-layer (0, 1, 2, 3, ...)into binary:
    
    
    # the user may want to drop constant columns: 
    if drop_const: 
        constant_columns = [col for col in data.columns if data[col].nunique() == 1]
        if verbose: print(f"WARNING: removing {len(constant_columns)} constant features.")
        data_bool = data_bool.drop(columns=constant_columns)

    
    # pdist() / linkage() will loose the accession information. So here we save a dict: 
    index_to_acc = {i: accession for i, accession in enumerate(data_bool.index)}
    # Calculate the linkage matrix using Ward clustering and Jaccard dissimilarity
    distances = pdist(data_bool, 'jaccard')
    linkage_matrix = linkage(distances, method='ward')
    
    
    # creates empty plots: 
    fig, axs = create_silhuette_frame(figsize)
    
    
    # get the vector of number of clusters to test:
    num_clusters_vector = np.arange(2, len(data_bool)-1, 1)
    if ctotest != None: num_clusters_vector = ctotest
    #print("Testing the following number of clusters:", num_clusters_vector)
    
    
    # Initialize lists to store silhouette scores and cluster assignments
    silhouette_avg_scores = []
    cluster_assignments = []
    
    
    # Iterate over a range of threshold values
    for num_clusters in num_clusters_vector:
        # Extract clusters based on the current threshold
        clusters = cut_tree(linkage_matrix, n_clusters=num_clusters)
        clusters = clusters.flatten()
        # 'clusters' is now a list of int, representing the cluster to which the i-element belongs to.
        # create a conversion dictionary: 
        acc_to_cluster = {index_to_acc[index]: clusters[index] for index in index_to_acc.keys()}
        
        
        # Calculate the silhouette score for the current set of clusters.
        # The Silhouette Score can be used for both K-means clustering and hierarchical clustering, 
        # as well as other clustering algorithms. It's a general-purpose metric for evaluating the 
        # quality of clusters, and it does not depend on the specific clustering algorithm being used.
        silhouette_avg = silhouette_score(data_bool, clusters)
        
        # Store the silhouette score and cluster assignments
        silhouette_avg_scores.append(silhouette_avg)
        cluster_assignments.append(clusters)

        
    # get the max average sillhouette (optimal value)
    max_value = max(silhouette_avg_scores)
    max_index = silhouette_avg_scores.index(max_value)
    opt_n_clusters = max_index + 2  # '+2' because num_clsuters starts from 2
    
        
    # Plot the average sihoutte (average on each datapoint).
    make_plot_ncluster_comparison(axs[0], num_clusters_vector, silhouette_avg_scores, opt_n_clusters, forcen, verbose)
    
        
    # Given the optimal number of clusters, visualizze the silhouette score for each data point. 
    if forcen != None: opt_n_clusters = forcen
    clusters = cut_tree(linkage_matrix, n_clusters=opt_n_clusters)
    clusters = clusters.flatten()
    acc_to_cluster = {index_to_acc[index]: clusters[index]+1 for index in index_to_acc.keys()}
    silhouette_avg = silhouette_score(data_bool, clusters)
    if verbose: print(f'Avg silhouette score when {opt_n_clusters} clusters:', silhouette_avg)
    silhouette_scores = silhouette_samples(data_bool, clusters)
    # Now 'silhouette_scores' is just a list of values. But the index correspond to a specific accession, that is
    # associated to a specific cluster. Thus, later we obtain the scores for a specific cluster 
    # simply with a 'silhouette_scores[clusters == i]'.
    
    
    # Show silhouette scores for each datapoint (given the opimal number of clusters)
    cluster_to_color = make_plot_silhouette_coeff(axs[2], opt_n_clusters, silhouette_scores, clusters)
    
    # Plot the dendrogram
    make_dendrogram(axs[4], linkage_matrix)
    
    # order the dataframe following the leaves of the tree:
    ord_data_bool = sort_by_leaves(data_bool, linkage_matrix, index_to_acc)
    
    # add colorbar for the dendrogram
    make_colorbar_clusters(axs[5], ord_data_bool, acc_to_cluster, cluster_to_color)
    
    # add colorbar for the species/niches
    make_colorbar_metadata(axs[7], ord_data_bool, derive_report, report_key,  excludekeys, key_to_color)
    
    # make legeneds
    make_legends(axs[9], derive_report, report_key,  excludekeys, cluster_to_color, None, anchor, key_to_color)
        
    # save to disk; bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
        
        
    return (fig, acc_to_cluster, cluster_to_color)



def heatmap_multilayer(
    tables, figsize = (10,5), drop_const=True, 
    derive_report=None, report_key='species', excludekeys=[], acc_to_cluster=None, cluster_to_color=None, 
    legend_ratio=0.7, label_ratio=0.02, outfile=None, verbose=False, anchor=[None, None, None], key_to_color=None,
    xlabels=False):
    """Create a phylo-metabolic dendrogram.
    
    Args:
        tables (pnd.DataFrame): feature tables with genome accessions are in columns and features are in rows. 
            Can also be a dictionary of feature tables (example: ``{'auxotrophies': aux_df, 'substrates': sub_df})``. 
            In this case, any number of tables (pandas.DataFrame) can be used. 
            For each table, genome accessions are in columns, features are in rows.
            Directly compatible tables are: `rpam.csv`, `cnps.csv`, and `aux.csv` (all produced by `gempipe derive`).
        figsize (int, int): width and height of the figure.
        drop_const (bool): if `True`, remove constant features.
        derive_report (pandas.DataFrame): report table for the generation of strain-specific GSMMs, made by `gempipe derive` in the output directory (`derive_strains.csv`). 
        report_key (str): name of the attribute (column) appearing in `derive_report`, to be compared to the metabolilc clusters.
            Usually it is 'species' or 'niche'.
        excludekeys (list): keys (iches/species) not to show in legend. 
            Bug: no more than 1 key is allowed.
        acc_to_cluster (dict):  genome-to-cluster associations produced by `silhouette_analysis()`.
        cluster_to_color (dict):  cluster-to-RGB color associations produced by `silhouette_analysis()`.
        legend_ratio (float): space reserved for the legend.
        label_ratio (float): space reserved for the Y-axis labels.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages
        anchor (list): list of tuples (X,Y) for customixing the position of legends. 
            ``None`` will leave default positioning.
        key_to_color (dict): dict mapping each category in `report_key` to a color in the format ([0:1],[0:1],[0:1]).
            ``None`` will leave default color and order in the legend. 
        xlabels (bool): if `True`, show x-axis labels (feature IDs).
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the phylometabolic tree and associated heatmap.
            - pnd.DataFrame: table representing the binary features contained in the heatmap.
    """
    
    
    def create_heatmap_frame(figsize):
    
        # create the subplots: 
        fig, axs = plt.subplots(
            nrows=1, ncols=8, 
            figsize=figsize, # global dimensions.
            gridspec_kw={'width_ratios': [0.3, 0.04, 0.02, 0.04,  0.02, 0.94, label_ratio, legend_ratio ]}) # suplots width proportions. 
        # adjust the space between subplots: 
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[2].axis('off')  # remove frame and axis
        axs[4].axis('off')  # remove frame and axis
        axs[6].axis('off')  # remove frame and axis

        return fig, axs
    
    
    
    def make_plot_heatmap_multilayer(ax, ord_data, xlabels ):
                
        ax.matshow(
            ord_data,  
            cmap='viridis',
            vmin=ord_data.min().min(), vmax=ord_data.max().max(), # define ranges for the colormap.
            aspect='auto', # fixed axes and aspect adjusted to fit data.
            interpolation='none') # no interp. performed on Agg-ps-pdf-svg backends.

        # set x labels
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        if xlabels:   # show x-axis labels (feats IDs)
            ax.get_xaxis().set_visible(True)
            ax.set_xticks(range(len(ord_data.columns)))
            ax.set_xticklabels(list(ord_data.columns))
            ax.xaxis.set_ticks_position('bottom')
        
 
    
    # START
    # format input tables: 
    data, dict_tables = merge_tables(tables)
    data_bool = data.astype(bool)   # convert multi-layer (0, 1, 2, 3, ...)into binary:

    
    # the user may want to drop constant columns: 
    if drop_const: 
        constant_columns = [col for col in data.columns if data[col].nunique() == 1]
        if verbose: print(f"WARNING: removing {len(constant_columns)} constant features.")
        data      = data.drop(columns=constant_columns)
        data_bool = data_bool.drop(columns=constant_columns)
    
    
    # pdist() / linkage() will loose the accession information. So here we save a dict: 
    index_to_acc = {i: accession for i, accession in enumerate(data.index)}
    
    
    # create a dendrogram based on the jaccard distancies (dissimilarities): 
    distances = pdist(data_bool, 'jaccard')
    linkage_matrix = linkage(distances, method='ward')
    
    
    # create the empty figure frame:
    fig, axs = create_heatmap_frame(figsize)

    # plot the dendrogram
    make_dendrogram(axs[0], linkage_matrix)

    # order the dataframe following the leaves of the tree:
    ord_data = sort_by_leaves(data, linkage_matrix, index_to_acc)  
    ord_data_bool = sort_by_leaves(data_bool, linkage_matrix, index_to_acc)  # only to return 
    
    # plot the heatmap:
    make_plot_heatmap_multilayer(axs[5], ord_data, xlabels)
    
    # add the cluster information (coming from the silhouette analysis);
    make_colorbar_clusters(axs[1], ord_data, acc_to_cluster, cluster_to_color)
    
    # colorbar for the species/niche
    make_colorbar_metadata(axs[3], ord_data, derive_report, report_key,  excludekeys, key_to_color)
    
    # make legends
    make_legends(axs[7], derive_report, report_key,  excludekeys, cluster_to_color, dict_tables, anchor, key_to_color)
    
    # save to disk; bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
                
        
    return fig, ord_data_bool



def discriminant_feat(binary_feats, acc_to_cluster, cluster_to_color, threshold=0.90):
    """Extract discriminant features from cluster of strains.
    
    Args:
        tables (pnd.DataFrame): binary features table such as the one produced by `heatmap_multilayer`
            (genomes in row, binary featuresin column).
        acc_to_cluster (dict): dictionary such as the one produced by `silhouette_analysis``
            (accessions as keys, cluster assignment as value).
        cluster_to_color (dict): dictionary such as the one produced by `silhouette_analysis``
            (clusters as keys, colors as value).
        threshold (float): features are shown if at least one cluster has relative frequency
            >= `threshold` and, at the same time, at least another cluster has relative frequency
            <= 1-`threshold`.
        
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the discriminative binary features.
    """
    
    def get_contingency(binary_feats, feat_id):
        contingency_table = pnd.crosstab(
            binary_feats[feat_id],
            binary_feats['y'],
            margins = False)
        
        # limit case: '0' or '1' is mising: 
        if 0 not in contingency_table.index:
            new_row = pnd.DataFrame([[0] * contingency_table.shape[1]], columns=contingency_table.columns, index=[0])
            contingency_table = pnd.concat([new_row, contingency_table])
        if 1 not in contingency_table.index:
            new_row = pnd.DataFrame([[0] * contingency_table.shape[1]], columns=contingency_table.columns, index=[1])
            contingency_table = pnd.concat([new_row, contingency_table])
        
        # the resulting pnd.DataFrame will be similar to (e.g. for the binary feat "[aux]his__L"):
        #         y            Cluster_1  Cluster_2  Cluster_3  Cluster_4  Cluster_5
        # [aux]his__L                                                       
        # 0                   12          3          6         12          0
        # 1                    0          0          0          0          3
        return contingency_table
    
    
    # START
    # convert to int:
    binary_feats = binary_feats.copy().astype(int)  # .copy() will defragment the dataframe.
    
    
    # add the classification column (usually caled 'y'):
    binary_feats['y'] = "Cluster_" + str(0)
    for accession, row in binary_feats.iterrows(): 
        binary_feats.loc[accession, 'y'] = "Cluster_" + str(acc_to_cluster[accession]) # str to avoid ambiguity
        
    
    # get dataframe of relative frequencies:
    df_relfreq = pnd.DataFrame(index=list(set(list(binary_feats.columns))-set(['y'])), columns=binary_feats['y'].unique())
    for feat_id in df_relfreq.index:
        cont = get_contingency(binary_feats, feat_id)
        for cluster in df_relfreq.columns:
            df_relfreq.loc[feat_id, cluster] = cont.loc[1, cluster] / (cont.loc[1, cluster] + cont.loc[0, cluster])
    df_relfreq  = df_relfreq.astype(float)
    
    
    # filter the dataframe:
    df_relfreq = df_relfreq[(df_relfreq >= threshold).any(axis=1)]
    df_relfreq = df_relfreq[(df_relfreq <= 1-threshold).any(axis=1)]
        
        
    # invert column order to match that of the heatmap
    df_relfreq = df_relfreq[reversed(df_relfreq.columns)]
    
    
    # sort features aalphabetically
    df_relfreq = df_relfreq.sort_index(ascending=False)
    
    # resort index in order to have similar rows (similar features) close together.
    index_to_featid = {i: feat_id for i, feat_id in enumerate(df_relfreq.index)}
    distances = pdist(df_relfreq, 'jaccard')
    linkage_matrix = linkage(distances, method='ward')
    df_relfreq = sort_by_leaves(df_relfreq, linkage_matrix, index_to_featid)  # only to return 
    
    
    # create the subplots: 
    fig, axs = plt.subplots(
        nrows=2, ncols=1, 
        figsize=(0.5 * len(df_relfreq.columns), 0.3 * len(df_relfreq)), # global dimensions.
        gridspec_kw={'width_ratios': [1], 'height_ratios': [1, 1*len(df_relfreq)]}) # suplots width proportions. 
    # adjust the space between subplots: 
    plt.subplots_adjust(wspace=0, hspace=0)
    axs[0].set_frame_on(False)  # remove squared border but not ticks
    axs[1].set_frame_on(False)  # remove squared border but not ticks
    
    
    # create matshow (clusters)
    if type(list(cluster_to_color.keys())[0]) == str: 
        cmap = LinearSegmentedColormap.from_list('', [cluster_to_color[cluster.replace('Cluster_', '')] for cluster in df_relfreq.columns])
    if type(list(cluster_to_color.keys())[0]) == int:   # add the 'eval'
        cmap = LinearSegmentedColormap.from_list('', [cluster_to_color[eval(cluster.replace('Cluster_', ''))] for cluster in df_relfreq.columns])
    df_clusters = pnd.DataFrame({cluster: [i] for i, cluster in enumerate(df_relfreq.columns)})
    axs[0].matshow(
        df_clusters,  
        cmap=cmap,
        vmin=df_clusters.min().min(), vmax=df_clusters.max().max(), # define ranges for the colormap.
        aspect='auto', # fixed axes and aspect adjusted to fit data.
        interpolation='none') # no interp. performed on Agg-ps-pdf-svg backends.
    
    
    # create matshow (heatmap)
    cmap = LinearSegmentedColormap.from_list('', ["#DDDDDD", "#8888DD"])
    axs[1].matshow(
        df_relfreq,  
        cmap=cmap,
        vmin=df_relfreq.min().min(), vmax=df_relfreq.max().max(), # define ranges for the colormap.
        aspect='auto', # fixed axes and aspect adjusted to fit data.
        interpolation='none') # no interp. performed on Agg-ps-pdf-svg backends.
    
    
    # add x/y axis labels (clusters):
    axs[0].set_xticks(range(len(df_relfreq.columns)))  # Position of x-ticks
    axs[0].set_xticklabels(df_relfreq.columns, ha='left', rotation=30)  # Replace with your desired labels
    axs[0].set_yticks([])   # remove x ticks
    
    
    # add x/y axis labels (heatmap):
    axs[1].set_xticks([])   # remove x ticks
    axs[1].set_yticks(range(len(df_relfreq.index)))  # Position of y-ticks
    axs[1].set_yticklabels(df_relfreq.index)
    
    
    # add annotations (rel frequencies):
    for i in range(len(df_relfreq.index)):
        for j in range(len(df_relfreq.columns)):
            freq_to_show = round(df_relfreq.iloc[i,j], 2)
            annot_color = 'white' if freq_to_show >= 0.60 else 'black'
            axs[1].text(j, i, f'{freq_to_show}', ha='center', va='center', color=annot_color)
    

    return fig, df_relfreq