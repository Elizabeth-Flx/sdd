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
from .utils import *
from .plots import *
import argparse

parser = argparse.ArgumentParser(description="Metrics evaluation for clustering methods.")

parser.add_argument("--input_anndata", type=str, required=True, help="Input anndata file directory")
parser.add_argument("--input_csv", type=str, required=True, help="Input csv files directory")
parser.add_argument("--output", type=str, required=True, help="Output figures directory")    

#parser.add_argument("--n_clusters", type=int, default=12, help="Number of clusters for clustering.")

args = parser.parse_args()

anndata_path = args.input_anndata
csv_path = args.input_csv
figures_path = args.output
print('figures_path:', figures_path)

#dirs = [
#    d for d in os.listdir(data_path)
#    if os.path.isdir(os.path.join(data_path, d)) and d != "figures"]
csv_dir = os.path.abspath(csv_path)
anndata_dir = os.path.abspath(anndata_path)
figures_dir = os.path.abspath(figures_path)

# make directory data_path and figures if it doesnt exist yet
#figures_path = os.path.join(data_path, "figures")
#print(f"Figures will be saved to: {figures_dir}")
#if not os.path.exists(figures_dir):
#    os.makedirs(figures_dir)

ARI_all = {}
Purity_all = {}
AMI_all = {}

ILISI_all = {}
PAS_all = {}

Silhouette_exp_all = {}
Silhouette_spat_all = {}

    

#print(f"Directory found: {data_path_dir}")

#data_path_dir = os.path.join(data_path, figures_dir)
#dirs_in_dir = os.listdir(data_path_dir)

dict_all_methods, groundtruth = load_csvs_to_arrays(csv_dir, os.listdir(csv_dir))

# adata as file in path that ends with .h5ad
for file in os.listdir(anndata_dir):
    if file.endswith(".h5ad"):
        h5ad_file = os.path.join(anndata_dir, file)
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
umap_path = os.path.join(figures_dir, f"umaps")
if not os.path.exists(umap_path):
    os.makedirs(umap_path)
    
umap_plot(adata, dict_all_methods, umap_path, groundtruth=groundtruth)

print(f"UMAP plots saved to {umap_path}")

sankey_path = os.path.join(figures_dir, f"sankey")
if not os.path.exists(sankey_path):
    os.makedirs(sankey_path)
    
fig_sankey = plot_sankey_diagram(dict_all_methods, ground_truth=groundtruth,
                                title="Sankey Diagram of Clustering Methods vs Ground Truth")

#fig_sankey.show()
fig_sankey.write_image(os.path.join(sankey_path, "sankey_diagram.png"))

print(f"Sankey diagram saved to {sankey_path}")
    
    
# accuracy metrics calculation
ARI_scores, Purity_scores, AMI_scores = calculate_accuracy_metrics(dict_all_methods, groundtruth=groundtruth)
ARI_all[figures_dir] = ARI_scores
Purity_all[figures_dir] = Purity_scores
AMI_all[figures_dir] = AMI_scores

print(f"Accuracy metrics calculated for directory {figures_dir}")

# spatial metrics calculation
ILISI_scores, PAS_scores = calculate_spatial_metrics(adata, dict_all_methods)
ILISI_all[figures_dir] = ILISI_scores
PAS_all[figures_dir] = PAS_scores

print(f"Spatial metrics calculated for directory {figures_dir}")

# silhouette metrics calculation
Silhouette_exp_scores, Silhouette_spat_scores = calculate_silhouette_metrics(adata, dict_all_methods)
Silhouette_exp_all[figures_dir] = Silhouette_exp_scores
Silhouette_spat_all[figures_dir] = Silhouette_spat_scores

print(f"Silhouette metrics calculated for directory {figures_dir}")
    
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

df_all_metrics.to_csv(os.path.join(figures_dir, "all_metrics_detailed.csv"))
df_all_metrics_mean.to_csv(os.path.join(figures_dir, "all_metrics_summary.csv"))

print(f"All metrics dataframes saved to {figures_dir}")

# here everything for plotting should be included in shiny as options and only then plotted
# plotting all
print(df_all_metrics_mean)
df_metrics_accuracy = df_all_metrics_mean[df_all_metrics_mean.columns.intersection(['Method', 'ARI', 'Purity', 'AMI'])]
print(df_metrics_accuracy)
df_metrics_spatial = df_all_metrics_mean[df_all_metrics_mean.columns.intersection(['Method', 'LISI', 'PAS'])]
print(df_metrics_spatial)
df_metrics_silhouette = df_all_metrics_mean[df_all_metrics_mean.columns.intersection(['Method', 'Silhouette_exp', 'Silhouette_spat'])]
boxplot_p = plot_metric_boxplots_gg(df_all_metrics, 'ARI', figures_dir)
lineplot_p = plot_metric_lineplot_datasets_gg(df_all_metrics, 'ARI', figures_dir)

mean_plots = plot_means_lineplot(df_all_metrics_mean, figures_dir, title="Mean Metric Scores Method")
mean_plots_accuracy = plot_means_lineplot(df_metrics_accuracy, figures_dir, title="Mean Metric Scores Accuracy")
mean_plots_spatial = plot_means_lineplot(df_metrics_spatial, figures_dir, title="Mean Metric Scores Spatial")
mean_plots_silhouette = plot_means_lineplot(df_metrics_silhouette, figures_dir, title="Mean Metric Scores Silhouette")

print(f"Plots generated and saved for directory {figures_dir}")



