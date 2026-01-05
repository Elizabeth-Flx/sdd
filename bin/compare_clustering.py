import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
from matplotlib_venn import venn2
import os
import argparse as ap

def get_cluster_labels(adata):
    for key in ["leiden", "louvain", "mclust", "cluster", "clusters"]:
        if key in adata.obs:
            return adata.obs[key].astype(str)
    # Fallback: use first categorical column
    for col in adata.obs.columns:
        if str(adata.obs[col].dtype) == "category" or adata.obs[col].dtype == object:
            return adata.obs[col].astype(str)
    raise ValueError("No cluster annotation found in AnnData.obs.")

def cluster_stats(y):
    sizes = pd.Series(y).value_counts().sort_index()
    proportions = sizes / sizes.sum()
    return sizes, proportions

def plot_cluster_sizes(sizes, out_path, name):
    plt.figure(figsize=(6,4))
    sizes.plot(kind='bar')
    plt.title(f"Cluster Size Distribution: {name}")
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{out_path}/cluster_size_{name}.png")
    plt.close()

def plot_umap(adata, y, out_path, name):
    if "X_umap" in adata.obsm:
        coords = adata.obsm["X_umap"]
        plt.figure(figsize=(7,6))
        plt.scatter(coords[:,0], coords[:,1], c=pd.Series(y).astype('category').cat.codes.to_numpy(), cmap='tab20', s=10)
        plt.title(f"UMAP: {name}")
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.tight_layout()
        plt.savefig(f"{out_path}/umap_{name}.png")
        plt.close()

def marker_genes(adata, y, top_n=10):
    markers = {}
    for cl in pd.Series(y).unique():
        mask = (y == cl)
        means = adata[mask].X.mean(axis=0) - adata[~mask].X.mean(axis=0)
        top = np.array(adata.var_names)[np.argsort(means)[-top_n:]]
        markers[cl] = set(top)
    return markers

def plot_marker_heatmap(adata, y, markers, out_path, name):
    all_markers = list(set.union(*markers.values()))
    sc.pl.heatmap(adata, all_markers, y.name) # pyright: ignore[reportCallIssue]
    plt.tight_layout()
    plt.savefig(f"{out_path}/marker_heatmap_{name}.png")
    plt.close()

def marker_overlap(markers1, markers2, out_path):
    # Flatten marker sets
    set1 = set.union(*markers1.values())
    set2 = set.union(*markers2.values())
    plt.figure(figsize=(5,5))
    venn2((set1, set2), set_labels=('Dataset1', 'Dataset2'))
    plt.title("Marker Gene Overlap")
    plt.tight_layout()
    plt.savefig(f"{out_path}/marker_gene_overlap.png")
    plt.close()
    # Jaccard index
    jaccard = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else np.nan
    with open(f"{out_path}/marker_gene_overlap.txt", "w") as f:
        f.write(f"Jaccard index: {jaccard}\n")
        f.write(f"Shared markers: {set1 & set2}\n")

def silhouette(adata, y, out_path, name):
    if "X_pca" in adata.obsm:
        sil = silhouette_score(adata.obsm['X_pca'], y)
    else:
        sil = np.nan
    with open(f"{out_path}/silhouette_{name}.txt", "w") as f:
        f.write(f"Silhouette score: {sil}\n")

def main():
    # Replace with argparse or manual paths

    parser = ap.ArgumentParser(description="Compare clustering results between two AnnData datasets.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing input AnnData files")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save comparison results")

    args = parser.parse_args()

    # Make output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    sp_file = os.path.join(args.input_dir, "c_sp_ad_graphst.h5ad")
    atac_file = os.path.join(args.input_dir, "sp_atac_graphst.h5ad")

    ad_sp = sc.read_h5ad(sp_file)
    ad_atac = sc.read_h5ad(atac_file)

    y_sp = get_cluster_labels(ad_sp)
    y_atac = get_cluster_labels(ad_atac)
    
    # Cluster stats
    sizes1, props1 = cluster_stats(y_sp)
    sizes2, props2 = cluster_stats(y_atac)
    plot_cluster_sizes(sizes1, args.output_dir, "sp")
    plot_cluster_sizes(sizes2, args.output_dir, "atac")
    # UMAP plots
    plot_umap(ad_sp, y_sp, args.output_dir, "sp")
    plot_umap(ad_atac, y_atac, args.output_dir, "atac")

    # Marker genes
    markers1 = marker_genes(ad_sp, y_sp, top_n=10)
    markers2 = marker_genes(ad_atac, y_atac, top_n=10)
    plot_marker_heatmap(ad_sp, y_sp, markers1, args.output_dir, "sp")
    plot_marker_heatmap(ad_atac, y_atac, markers2, args.output_dir, "atac")
    marker_overlap(markers1, markers2, args.output_dir)

    # Silhouette scores
    silhouette(ad_sp, y_sp, args.output_dir, "sp")
    silhouette(ad_atac, y_atac, args.output_dir, "atac")

    # Save summary stats
    with open(f"{args.output_dir}/cluster_stats.txt", "w") as f:
        f.write("Dataset1 cluster sizes:\n")
        f.write(str(sizes1) + "\n")
        f.write("Dataset1 cluster proportions:\n")
        f.write(str(props1) + "\n")
        f.write("Dataset2 cluster sizes:\n")
        f.write(str(sizes2) + "\n")
        f.write("Dataset2 cluster proportions:\n")
        f.write(str(props2) + "\n")

if __name__ == "__main__":
    main()