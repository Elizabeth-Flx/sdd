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
    # Prefer mclust if available
    if "mclust" in adata.obs:
        return adata.obs["mclust"].astype(str)
    # Fallbacks...
    for key in ["leiden", "louvain", "cluster", "clusters", "domain"]:
        if key in adata.obs:
            return adata.obs[key].astype(str)
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
    # Convert index to int for sorting, then back to string for labeling
    try:
        cluster_numbers = [int(str(c)) for c in sizes.index]
    except Exception:
        # If conversion fails, fallback to original order
        cluster_numbers = list(range(1, len(sizes.index)+1))
    order = np.argsort(cluster_numbers)
    ordered_sizes = sizes.iloc[order]
    # Rename to cl_<clusternumber>
    new_labels = [f"cl_{int(str(c))}" for c in ordered_sizes.index]
    plt.figure(figsize=(6,4))
    bars = plt.bar(new_labels, ordered_sizes.values, edgecolor='black', linewidth=1, width=1.0)
    plt.title(f"Cluster Size Distribution: {name}")
    plt.xlabel('Cluster')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{out_path}/cluster_size_{name}.png")
    plt.close()

def plot_umap_from_obs(adata, y, out_path, name, umap_x="HVG_UMAP_spatial", umap_y=None):
    if umap_x in adata.obs:
        if umap_y is None:
            # Try to infer the y column
            if umap_x.endswith("_x"):
                umap_y = umap_x[:-2] + "_y"
            else:
                print(f"Cannot infer UMAP y column for {umap_x}")
                return
        x = adata.obs[umap_x]
        y_coord = adata.obs[umap_y]
        cluster_codes = pd.Series(y).astype('category').cat.codes.to_numpy()
        plt.figure(figsize=(7,6))
        plt.scatter(x, y_coord, c=cluster_codes, cmap='tab20', s=10)
        plt.title(f"UMAP: {name}")
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.tight_layout()
        plt.savefig(f"{out_path}/umap_{name}.png")
        plt.close()

def marker_genes(adata, y, top_n=10):
    markers = {}
    # print type of y
    print("Type of y:", type(y))
    y_series = pd.Series(y)
    for cl in y_series.unique():
        mask = (y_series == cl)
        # Ensure means is a 1D array
        means = adata[mask].X.mean(axis=0) - adata[~mask].X.mean(axis=0)
        if hasattr(means, 'A1'):  # handle sparse matrix
            means = means.A1
        means = np.array(means).flatten()
        var_names = np.array(adata.var_names).flatten()
        # Select top N gene indices
        top_indices = np.argsort(means)[-top_n:]
        top_genes = var_names[top_indices]
        top_genes = top_genes.tolist() if isinstance(top_genes, np.ndarray) else list(top_genes)
        print(f"Cluster {cl} top markers: {len(top_genes)} genes: {top_genes}")
        markers[cl] = set(top_genes)
    return markers


def plot_marker_heatmap(adata, y, markers, out_path, name):
    all_markers = list(set.union(*markers.values()))
    print(f"Total unique markers for {name}: {len(all_markers)}")
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



def main():
    # Replace with argparse or manual paths

    parser = ap.ArgumentParser(description="Compare clustering results between two AnnData datasets.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing input AnnData files")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save comparison results")

    args = parser.parse_args()

    # Make output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    dlpfc_file = os.path.join(args.input_dir, "sp_dlpfc_graphst.h5ad")
    atac_file = os.path.join(args.input_dir, "sp_atac_graphst.h5ad")
    sp_file = os.path.join(args.input_dir, "c_sp_ad_graphst.h5ad")

    ad_dlpfc = sc.read_h5ad(dlpfc_file)
    ad_atac = sc.read_h5ad(atac_file)
    ad_sp = sc.read_h5ad(sp_file)

    y_dlpfc = get_cluster_labels(ad_dlpfc)
    y_atac = get_cluster_labels(ad_atac)
    y_sp = get_cluster_labels(ad_sp)
    
    # Cluster stats
    sizes1, props1 = cluster_stats(y_dlpfc)
    sizes2, props2 = cluster_stats(y_atac)
    sizes3, props3 = cluster_stats(y_sp)
    plot_cluster_sizes(sizes1, args.output_dir, "dlpfc")
    plot_cluster_sizes(sizes2, args.output_dir, "atac")
    plot_cluster_sizes(sizes3, args.output_dir, "sp")
    # UMAP plots
    plot_umap_from_obs(ad_dlpfc, y_dlpfc, args.output_dir, "dlpfc", umap_x="HVG_UMAP_spatial_x", umap_y="HVG_UMAP_spatial_y")
    plot_umap_from_obs(ad_atac, y_atac, args.output_dir, "atac", umap_x="HVG_UMAP_spatial_x", umap_y="HVG_UMAP_spatial_y")
    plot_umap_from_obs(ad_sp, y_sp, args.output_dir, "sp", umap_x="HVG_UMAP_spatial_x", umap_y="HVG_UMAP_spatial_y")

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
    marker_overlap(markers1, markers2, args.output_dir)

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
        f.write("Dataset3 cluster sizes:\n")
        f.write(str(sizes3) + "\n")
        f.write("Dataset3 cluster proportions:\n")
        f.write(str(props3) + "\n")

if __name__ == "__main__":
    main()