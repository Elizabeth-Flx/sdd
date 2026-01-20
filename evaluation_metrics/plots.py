import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import *
from plotnine import *
import scanpy as sc
import anndata as ad
from matplotlib.path import Path
import os
from scib.metrics import *
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
import plotly.graph_objects as go
import numpy as np
import pandas as pd

def umap_plot(adata, dict_all_methods, umap_path, groundtruth=None):
    
    sc.settings.figdir = Path(umap_path)

    if not os.path.exists(umap_path):
        os.makedirs(umap_path)

    if groundtruth is not None:
            adata.obs['groundtruth'] = pd.Categorical(groundtruth)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.pl.umap(adata, color='groundtruth', title='groundtruth', save = "_groundtruth.png")
            # save umap plot
            #plt.savefig(os.path.join(umap_path, "umap_groundtruth.png"))
            #plt.close() 
            
    for method_name, values in dict_all_methods.items():
        adata.obs[method_name] = pd.Categorical(values)
        if 'X_umap' not in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        sc.pl.umap(adata, color=method_name, title=method_name, show=False, save=f"_{method_name}.png")
        # save umap plot
        #plt.savefig(os.path.join(umap_path, f"umap_{method_name}.png"))
        #plt.close()



# sankey plot function

def plot_sankey_diagram(methods_dict, ground_truth=None, 
                        method_colors=None, title="Sankey Diagram of Clustering Methods"):
    """
    Create Sankey diagram showing relationships between clustering methods
    
    Parameters:
    -----------
    methods_dict : dict
        Dictionary with method names as keys and cluster labels as values
        Example: {'EnSDD': ensdd_labels, 'STAGATE': stagate_labels, 'SpaGCN': spagcn_labels}
    ground_truth : array-like, optional
        Ground truth labels (expert annotation)
    method_colors : dict, optional
        Dictionary mapping method names to colors
    title : str
        Title of the plot
    """
    
    # Add ground truth if provided
    if ground_truth is not None:
        methods_dict = {'Ground Truth': ground_truth, **methods_dict}
    
    # Prepare data
    method_names = list(methods_dict.keys())
    n_methods = len(method_names)
    n_spots = len(list(methods_dict.values())[0])
    
    # Create all node labels: "Method_Cluster#"
    node_labels = []
    node_methods = []  # Track which method each node belongs to
    node_indices = {}  # Map from (method, cluster) to node index
    
    current_index = 0
    for method_idx, method_name in enumerate(method_names):
        labels = methods_dict[method_name]
        unique_clusters = np.unique(labels)
        
        for cluster in unique_clusters:
            node_label = f"{method_name}: {cluster}"
            node_labels.append(node_label)
            node_methods.append(method_name)
            node_indices[(method_name, cluster)] = current_index
            current_index += 1
    
    n_nodes = len(node_labels)
    
    # Create links between consecutive methods
    source = []
    target = []
    value = []
    link_colors = []
    
    # Define default colors for methods
    if method_colors is None:
        # Use Plotly's categorical colors (without transparency)
        default_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]
        method_colors = {method: default_colors[i % len(default_colors)] 
                        for i, method in enumerate(method_names)}
    
    # Create links between consecutive pairs of methods
    for i in range(n_methods - 1):
        method_a = method_names[i]
        method_b = method_names[i + 1]
        
        labels_a = methods_dict[method_a]
        labels_b = methods_dict[method_b]
        
        # Create contingency matrix
        unique_a = np.unique(labels_a)
        unique_b = np.unique(labels_b)
        
        for idx_a, cluster_a in enumerate(unique_a):
            for idx_b, cluster_b in enumerate(unique_b):
                # Count spots that are in cluster_a in method_a and cluster_b in method_b
                count = np.sum((labels_a == cluster_a) & (labels_b == cluster_b))
                
                if count > 0:  # Only add if there are spots
                    source_idx = node_indices[(method_a, cluster_a)]
                    target_idx = node_indices[(method_b, cluster_b)]
                    
                    source.append(source_idx)
                    target.append(target_idx)
                    value.append(count)
                    
                    # Color based on source method - WITHOUT TRANSPARENCY
                    link_colors.append(method_colors[method_a])
    
    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels,
            color=[method_colors[method] for method in node_methods],
            hovertemplate='%{label}<br>Total spots: %{value}<extra></extra>'
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=link_colors,
            hovertemplate='From %{source.label}<br>To %{target.label}<br>Spots: %{value}<extra></extra>'
        )
    )])
    
    # Update layout
    fig.update_layout(
        #title_text=title,
        title_font_size=20,
        font_size=12,
        height=600,
        annotations=[
            dict(
                x=0.5,
                y=1.1,
                xref="paper",
                yref="paper",
                text=title,
                showarrow=False,
                font=dict(size=16)
            )
        ]
    )
    
    return fig

# Boxplots for one metric for all methods and datasets
def plot_metric_boxplots_gg(metric_df, metric_name, save_path):
    
    # Create boxplot using plotnine (ggplot2-like)
    # also include all points as jittered points
    p = (ggplot(metric_df, aes(x='Method', y=metric_name, fill = 'Method')) +
         geom_boxplot() +
         geom_jitter(width=0.2, alpha=0.5, color='blue') +
         ggtitle(f'Boxplot of {metric_name} Scores by Method') +
         xlab('Method') +
         ylab(metric_name) +
         theme(axis_text_x=element_text(rotation=45, hjust=1), legend_position='none'))
    
    
    # Save plot
    p.save(os.path.join(save_path, f'{metric_name}_boxplot.png'))
    return p

# 

# Lineplots 
def plot_metric_lineplot_datasets_gg(metric_df, metric_name, save_path):
    p = (ggplot(metric_df, aes(x='Dataset', y=metric_name, color=metric_name, group='Method')) +
         geom_line(size=1) +
         geom_point(size=2) +
         ggtitle(f'Line Plot of {metric_name} Scores by Method across Datasets') +
         xlab('Dataset') +
         ylab(metric_name) +
         theme(axis_text_x=element_text(rotation=45, hjust=1)))
    
    #p.show()
    p.save(os.path.join(save_path, f'{metric_name}_lineplot_dataset.png'))
    return p

# mean lineplots of all metrics per datasets
def plot_means_lineplot(mean_df, title, save_path = None):
    # long format: Method | Metric | Mean_Score
    df_long = mean_df.melt(id_vars=["Method"], var_name="Metric", value_name="Mean_Score")

    p = (
        ggplot(df_long, aes(x="Method", y="Mean_Score", color="Metric", group="Metric"))
        + geom_line(size=1)
        + geom_point(size=2)
        + ggtitle(title)
        + xlab("Method")
        + ylab("Mean score")
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
    )

    # show + save
    #print(p)
    if save_path is not None:
        p.save(os.path.join(save_path, f"{title}.png"), dpi=300)

    return p
    