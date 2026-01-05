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
import sys
from scib.metrics import *
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
from utils import *
from plots import *

#data_path = "/nfs/home/students/s.reif/SysMed/SysBioSDD/evaluation_metrics/examples/"
data_path = sys.argv[1]

dirs = [
    d for d in os.listdir(data_path)
    if os.path.isdir(os.path.join(data_path, d)) and d != "figures"]

# make directory data_path and figures if it doesnt exist yet
figures_path = os.path.join(data_path, "figures")
print(f"Figures will be saved to: {figures_path}")
if not os.path.exists(figures_path):
    os.makedirs(figures_path)

ARI_all = {}
Purity_all = {}
AMI_all = {}

ILISI_all = {}
PAS_all = {}

Silhouette_exp_all = {}
Silhouette_spat_all = {}

    
for dir in dirs:
    print(f"Directory found: {dir}")
    
    data_path_dir = os.path.join(data_path, dir)
    dirs_in_dir = os.listdir(data_path_dir)
    
    dict_all_methods, groundtruth = load_csvs_to_arrays(data_path_dir, dirs_in_dir)
    
    # adata as file in path that ends with .h5ad
    for file in os.listdir(data_path_dir):
        if file.endswith(".h5ad"):
            h5ad_file = os.path.join(data_path_dir, file)
            print(f"Loading AnnData from {h5ad_file}")
            adata = sc.read_h5ad(h5ad_file)
            break
    else:
    # want error here
        raise FileNotFoundError("No .h5ad file found in the directory.")
    
    # checking the inputs: same size of labels and adata observations
    
    n = adata.n_obs
    for k, v in dict_all_methods.items():
        if len(v) != n:
            raise ValueError(f"{k} has {len(v)} labels, but adata has {n} observations")

    if groundtruth is not None and len(groundtruth) != n:
        raise ValueError(f"groundtruth has {len(groundtruth)} labels, but adata has {n} observations")

    # load all labels of each method in adata.obs
    
    for method, labels in dict_all_methods.items():
        adata.obs[method] = pd.Categorical(labels)

    if groundtruth is not None:
        adata.obs["groundtruth"] = pd.Categorical(groundtruth)
        
    # umaps and sankey plots
    umap_path = os.path.join(figures_path, f"umaps_{dir}")
    if not os.path.exists(umap_path):
        os.makedirs(umap_path)
        
    umap_plot(adata, dict_all_methods, umap_path, groundtruth=groundtruth)
    
    sankey_path = os.path.join(figures_path, f"sankey_{dir}")
    if not os.path.exists(sankey_path):
        os.makedirs(sankey_path)
        
    fig_sankey = plot_sankey_diagram(dict_all_methods, ground_truth=groundtruth,
                                 title="Sankey Diagram of Clustering Methods vs Ground Truth")

    #fig_sankey.show()
    fig_sankey.write_image(os.path.join(sankey_path, "sankey_diagram.png"))
        
        
    # accuracy metrics calculation
    ARI_scores, Purity_scores, AMI_scores = calculate_accuracy_metrics(dict_all_methods, groundtruth=groundtruth)
    ARI_all[dir] = ARI_scores
    Purity_all[dir] = Purity_scores
    AMI_all[dir] = AMI_scores
    
    # spatial metrics calculation
    ILISI_scores, PAS_scores = calculate_spatial_metrics(adata, dict_all_methods)
    ILISI_all[dir] = ILISI_scores
    PAS_all[dir] = PAS_scores
    
    # silhouette metrics calculation
    Silhouette_exp_scores, Silhouette_spat_scores = calculate_silhouette_metrics(adata, dict_all_methods)
    Silhouette_exp_all[dir] = Silhouette_exp_scores
    Silhouette_spat_all[dir] = Silhouette_spat_scores
    
# generate dataframe for all metrics

metric_dicts = {
    "ARI": ARI_all,
    "Purity": Purity_all,
    "AMI": AMI_all,
    "LISI": ILISI_all,
    "PAS": PAS_all,
    "Silhouette_exp": Silhouette_exp_all,
    "Silhouette_spat": Silhouette_spat_all
}

df_all_metrics, df_all_metrics_mean = generate_metrics_df_from_dicts(metric_dicts)

df_all_metrics.to_csv(os.path.join(figures_path, "all_metrics_detailed.csv"))
df_all_metrics_mean.to_csv(os.path.join(figures_path, "all_metrics_summary.csv"))

# here everything for plotting should be included in shiny as options and only then plotted

boxplot_p = plot_metric_boxplots_gg(df_all_metrics, 'ARI', figures_path)
lineplot_p = plot_metric_lineplot_datasets_gg(df_all_metrics, 'ARI', figures_path)

mean_plots = plot_means_lineplot(df_all_metrics_mean, figures_path)




        
    