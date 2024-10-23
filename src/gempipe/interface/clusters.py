import pandas as pnd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram, leaves_list
from sklearn.metrics import silhouette_score, silhouette_samples
import numpy as np    


from .clusters_utils import merge_tables, make_dendrogram, make_colorbar_clusters, make_colorbar_metadata, make_legends

    
    
def silhouette_analysis(tables, figsize = (10,5), ctotest=None, forcen=None, derive_report=None, report_key='species', legend_ratio=0.4, outfile=None, verbose=False):
    """Perform a silhuette analysis to detect the optimal number of clusters. 
    
    Args:
        tables (pnd.DataFrame): feture tables with genome accessions are in columns and features are in rows. 
            Can also be a dictionary of feature tables (example: ``{'auxotrophies': aux_df, 'substrates': sub_df})``. 
            In this case, any number of tables (pandas.DataFrame) can be used. 
            For each table, genome accessions are in columns, features are in rows.
            Directly compatible tables are: `rpam.csv`, `cnps.csv`, and `aux.csv` (all produced by `gempipe derive`).
        figsize (int, int): width and height of the figure.
        ctotest (list): number of clusters to test (example: ``[5,7,10]`` to test five, seven and ten clusters).
            If `None`, all the combinations from 2 to the number of accessions -1 will be used.
        forcen (int): force the number of cluster, otherwise the optimal number will picked up according to the sihouette value. 
        derive_report (pandas.DataFrame): report table for the generation of strain-specific GSMMs, made by `gempipe derive` in the output directory (`derive_strains.csv`). 
        report_key (str): name of the attribute (column) appearing in `derive_report`, to be compared to the metabolilc clusters.
            Usually it is 'species' or 'niche'.
        legend_ratio (float): space reserved for the legend.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages
    
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



    def make_plot_1_silhouette(ax, num_clusters_vector, silhouette_avg_scores, opt_n_clusters, forcen, verbose):

        # Plot the silhouette scores against the number of clusters (threshold values)
        ax.plot(num_clusters_vector, silhouette_avg_scores, marker='o')
        ax.set_xlabel('N clusters')
        ax.set_ylabel('Average Silhouette Score')
        ax.grid(True)

        if verbose: print(f"Optimal number of clusters: {opt_n_clusters}")
        ax.axvline(x=opt_n_clusters if forcen==None else forcen, color='red', linestyle='--')



    def make_plot_2_silhouette(ax, opt_n_clusters, silhouette_scores, clusters):

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
            ax.text(0, (y_lower + y_upper -1)/2, f'Cluster_{i+1}')
            y_lower = y_upper + -1  # no space between clusters

        ax.set_xlabel('Silhouette Coefficient')
        ax.set_ylabel('')
        ax.set_yticks([]) 
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_title('Silhouette Plot for {} Clusters'.format(opt_n_clusters))
        ax.set_facecolor('#f8f8f8')  # Light gray background

        return cluster_to_color



    # START:
    # format input tables:
    data, dict_tables = merge_tables(tables)
    data = data.astype(bool)   # convert multi-layer (0, 1, 2, 3, ...)into binary:
    
    
    # pdist() / linkage() will loose the accession information. So here we save a dict: 
    index_to_acc = {i: accession for i, accession in enumerate(data.index)}
    # Calculate the linkage matrix using Ward clustering and Jaccard dissimilarity
    distances = pdist(data, 'jaccard')
    linkage_matrix = linkage(distances, method='ward')
    
    
    # creates empty plots: 
    fig, axs = create_silhuette_frame(figsize)
    
    
    # get the vector of number of clusters to test:
    num_clusters_vector = np.arange(2, len(data)-1, 1)
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
        silhouette_avg = silhouette_score(data, clusters)
        
        # Store the silhouette score and cluster assignments
        silhouette_avg_scores.append(silhouette_avg)
        cluster_assignments.append(clusters)

        
    # get the max average sillhouette (optimal value)
    max_value = max(silhouette_avg_scores)
    max_index = silhouette_avg_scores.index(max_value)
    opt_n_clusters = max_index + 2  # '+2' because num_clsuters starts from 2
    
        
    # FIRST PLOT. Plot the average sihoutte (average on each datapoint).
    make_plot_1_silhouette(axs[0], num_clusters_vector, silhouette_avg_scores, opt_n_clusters, forcen, verbose)
    
        
    # Given the optimal number of clusters, visualizze the silhouette score for each data point. 
    if forcen != None: opt_n_clusters = forcen
    clusters = cut_tree(linkage_matrix, n_clusters=opt_n_clusters)
    clusters = clusters.flatten()
    acc_to_cluster = {index_to_acc[index]: clusters[index]+1 for index in index_to_acc.keys()}
    silhouette_avg = silhouette_score(data, clusters)
    if verbose: print(f'Avg silhouette score when {opt_n_clusters} clusters:', silhouette_avg)
    silhouette_scores = silhouette_samples(data, clusters)
    # Now 'silhouette_scores' is just a list of values. But the index correspond to a specific accession, that is
    # associated to a specific cluster. Thus, later we obtain the scores for a specific cluster 
    # simply with a 'silhouette_scores[clusters == i]'.
    

    # SECOND PLOT. Show silhouette scores for each datapoint (given the opimal number of clusters)
    cluster_to_color = make_plot_2_silhouette(axs[2], opt_n_clusters, silhouette_scores, clusters)
    
    # THIRD PLOT. Plot the dendrogram
    make_dendrogram(axs[4], linkage_matrix)
    
    # FOURTH PLOT: add colorbar for the dendrogram
    make_colorbar_clusters(axs[5], index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix)
    
    # FIFTH PLOT: add colorbar for the species/niches
    make_colorbar_metadata(axs[7], derive_report, report_key, index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix)
    
    # make legeneds
    make_legends(axs[9], derive_report, report_key, cluster_to_color, dict_tables=None)
        
    # save to disk; bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
        
        
    return (fig, acc_to_cluster, cluster_to_color)



def phylomet_dendro(tables, figsize = (10,5), drop_const=True, derive_report=None, report_key='species', acc_to_cluster=None, cluster_to_color=None, legend_ratio=0.4, outfile=None, verbose=False):
    """Create a phylo-metabolic dendrogram.
    
    Args:
        tables (pnd.DataFrame): feture tables with genome accessions are in columns and features are in rows. 
            Can also be a dictionary of feature tables (example: ``{'auxotrophies': aux_df, 'substrates': sub_df})``. 
            In this case, any number of tables (pandas.DataFrame) can be used. 
            For each table, genome accessions are in columns, features are in rows.
            Directly compatible tables are: `rpam.csv`, `cnps.csv`, and `aux.csv` (all produced by `gempipe derive`).
        figsize (int, int): width and height of the figure.
        drop_const (bool): if `True`, remove constant features.
        derive_report (pandas.DataFrame): report table for the generation of strain-specific GSMMs, made by `gempipe derive` in the output directory (`derive_strains.csv`). 
        report_key (str): name of the attribute (column) appearing in `derive_report`, to be compared to the metabolilc clusters.
            Usually it is 'species' or 'niche'.
        acc_to_cluster (dict):  genome-to-cluster associations produced by `silhouette_analysis()`.
        cluster_to_color (dict):  cluster-to-RGB color associations produced by `silhouette_analysis()`.
        legend_ratio (float): space reserved for the legend.
        outfile (str): filepath to be used to save the image. If `None` it will not be saved.
        verbose (bool): if `True`, print more log messages
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: figure representing the sinhouette analysis.
    """
    
    
    def create_dendro_frame(figsize):
    
        # create the subplots: 
        fig, axs = plt.subplots(
            nrows=1, ncols=8, 
            figsize=figsize, # global dimensions.
            gridspec_kw={'width_ratios': [0.3, 0.04, 0.02, 0.94, 0.02, 0.04,  0.02, legend_ratio ]}) # suplots width proportions. 
        # adjust the space between subplots: 
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[2].axis('off')  # remove frame and axis
        axs[4].axis('off')  # remove frame and axis
        axs[6].axis('off')  # remove frame and axis

        return fig, axs
    
    
    
    def make_plot_1_dendro(ax, ord_data):
        
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
    fig, axs = create_dendro_frame(figsize)
    

    # FIRST PLOT: plot the dendrogram
    make_dendrogram(axs[0], linkage_matrix)
    
    # SECOND PLOT: add the cluster information (coming from the silhouette analysis);
    make_colorbar_clusters(axs[1], index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix)
 
    
    # How to get the leaves order: 
    ord_leaves = leaves_list(linkage_matrix)
    ord_leaves = np.flip(ord_leaves)  # because leaves are returned in the inverse sense.
    ord_leaves = [index_to_acc[i] for i in ord_leaves]  # convert index as number to index as accession
    ord_data = data.loc[ord_leaves, :]  # # transposed and reordered dataframe.
     
    
    # THIRD PLOT: plot the heatmap:
    make_plot_1_dendro(axs[3], ord_data)
    
    # FOURTH PLOT: colorbar for the species/niche
    make_colorbar_metadata(axs[5], derive_report, report_key, index_to_acc, acc_to_cluster, cluster_to_color, linkage_matrix)
    
    # make legends
    make_legends(axs[7], derive_report, report_key, cluster_to_color, dict_tables)
    
    # save to disk; bbox_inches='tight' removes white spaces around the figure. 
    if outfile != None:
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
        
        
    return fig