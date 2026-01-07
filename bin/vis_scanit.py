import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse as ap

from scipy import cluster
from sklearn import naive_bayes

def get_cluster_labels(adata):
    if "scanit_louvain" in adata.obs:
        return adata.obs["scanit_louvain"].astype(str)
    raise ValueError("No cluster annotation found in AnnData.obs.")


def plot_spatial(adata, y, out_path, name):
    if "array_col" in adata.obs and "array_row" in adata.obs:
        x = adata.obs["array_col"]
        y_coord = adata.obs["array_row"]
    elif "centroid_x" in adata.obs and "centroid_y" in adata.obs:
        x = adata.obs["centroid_x"]
        y_coord = adata.obs["centroid_y"]
    else:
        print(f"No spatial coordinates found for {name}.")
        return
    cluster_codes = pd.Series(y).astype('category').cat.codes.to_numpy()
    plt.figure(figsize=(7,6))
    plt.scatter(x, y_coord, c=cluster_codes, cmap='tab20', s=10)
    plt.title(f"Spatial Clusters: {name}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(f"{out_path}/spatial_{name}.png")
    plt.close()

def cluster_stats(y):
    sizes = pd.Series(y).value_counts().sort_index()
    proportions = sizes / sizes.sum()
    return sizes, proportions


def plot_cluster_sizes(sizes, out_path, name):
    try:
        cluster_numbers = [int(str(c)) for c in sizes.index]
    except Exception:
        cluster_numbers = list(range(1, len(sizes.index)+1))
    order = np.argsort(cluster_numbers)
    ordered_sizes = sizes.iloc[order]
    new_labels = [f"cl_{int(str(c))}" for c in ordered_sizes.index]
    plt.figure(figsize=(6,4))
    bars = plt.bar(new_labels, ordered_sizes.values, edgecolor='black', linewidth=1, width=1.0)
    plt.title(f"Cluster Size Distribution: {name}")
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{out_path}/cluster_size_{name}.png")
    plt.close()

def plot_umap_from_obs(adata, out_path, name):
    out_file = f"{out_path}/umap_{name}.png"
    print(f"Saving UMAP plot to {out_file}")
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)  # pyright: ignore[reportCallIssue]
    figs = sc.pl.umap(adata, color='scanit_louvain', show=False, return_fig=True)
    # sc.pl.umap with return_fig=True returns a Figure or list of Figures
    if isinstance(figs, list):
        fig = figs[0]
    else:
        fig = figs
    fig.savefig(out_file)  # type: ignore
    plt.close()

def marker_genes(adata, y, top_n=10):
    markers = {}
    y_series = pd.Series(y)
    for cl in y_series.unique():
        mask = (y_series == cl)
        means = adata[mask].X.mean(axis=0) - adata[~mask].X.mean(axis=0)
        if hasattr(means, 'A1'):
            means = means.A1
        means = np.array(means).flatten()
        var_names = np.array(adata.var_names).flatten()
        top_indices = np.argsort(means)[-top_n:]
        top_genes = var_names[top_indices]
        top_genes = top_genes.tolist() if isinstance(top_genes, np.ndarray) else list(top_genes)
        markers[cl] = set(top_genes)
    return markers


def plot_marker_heatmap(adata, y, markers, out_path, name):
    all_markers = list(set.union(*markers.values()))
    sc.pl.heatmap(adata, all_markers, y.name) # pyright: ignore[reportCallIssue]
    plt.savefig(f"{out_path}/marker_heatmap_{name}.png")
    plt.close()


def main():
    parser = ap.ArgumentParser(description="Compare clustering results between two AnnData datasets.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing input AnnData files")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save comparison results")

    args = parser.parse_args()

    cluster_key = "scanit_louvain"

    # Make output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Use tool name in filenames
    dlpfc_file = os.path.join(args.input_dir, "sp_dlpfc_scanit.h5ad")
    atac_file = os.path.join(args.input_dir, "sp_atac_scanit.h5ad")
    sp_file = os.path.join(args.input_dir, "sp_scanit.h5ad")

    ad_dlpfc = sc.read_h5ad(dlpfc_file)
    ad_atac = sc.read_h5ad(atac_file)
    ad_sp = sc.read_h5ad(sp_file)
    print(ad_sp.obs.head())
    print(ad_sp.obs.columns)

    y_dlpfc = get_cluster_labels(ad_dlpfc)
    y_atac = get_cluster_labels(ad_atac)
    y_sp = get_cluster_labels(ad_sp)
    
    # UMAP plots
    plot_umap_from_obs(ad_dlpfc, args.output_dir, "dlpfc")
    plot_umap_from_obs(ad_atac, args.output_dir, "atac")
    plot_umap_from_obs(ad_sp, args.output_dir, "sp")

    # Spatial plots
    plot_spatial(ad_dlpfc, y_dlpfc, args.output_dir, "dlpfc")
    plot_spatial(ad_atac, y_atac, args.output_dir, "atac")
    plot_spatial(ad_sp, y_sp, args.output_dir, "sp")

    # Marker genes
    markers1 = marker_genes(ad_dlpfc, y_dlpfc, top_n=10)
    markers2 = marker_genes(ad_atac, y_atac, top_n=10)
    markers3 = marker_genes(ad_sp, y_sp, top_n=10)
    plot_marker_heatmap(ad_dlpfc, y_dlpfc, markers1, args.output_dir, "dlpfc")
    plot_marker_heatmap(ad_atac, y_atac, markers2, args.output_dir, "atac")
    plot_marker_heatmap(ad_sp, y_sp, markers3, args.output_dir, "sp")

if __name__ == "__main__":
    main()
