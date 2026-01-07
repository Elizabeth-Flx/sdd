from matplotlib.pylab import f
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse as ap


def get_cluster_labels(adata):
    # Prefer mclust if available
    if "stlearn_louvain" in adata.obs:
        return adata.obs["stlearn_louvain"].astype(str)
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


def plot_umap_from_obs(adata, color_col, out_path, name):
    out_file = f"{out_path}/umap_{name}.png"
    print(f"Saving UMAP plot to {out_file}")
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)  # pyright: ignore[reportCallIssue]
    figs = sc.pl.umap(adata, color=color_col, show=False, return_fig=True)
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


def display_kmeans_clustering(adata, out_path, name) -> None:
    if "kmeans" not in adata.obs:
        print(f"No kmeans clustering found for {name}.")
        return
    y_kmeans = adata.obs["kmeans"].astype(str)
    plot_spatial(adata, y_kmeans, out_path, f"{name}_kmeans")
    plot_umap_from_obs(adata, "kmeans", out_path, f"{name}_kmeans")
    sizes, _ = cluster_stats(y_kmeans)
    plot_cluster_sizes(sizes, out_path, f"{name}_kmeans")
    # markers = marker_genes(adata, y_kmeans, top_n=10)
    # plot_marker_heatmap(adata, y_kmeans, markers, out_path, f"{name}_kmeans")


def main():
    parser = ap.ArgumentParser(description="Visualize clustering and marker genes for stlearn AnnData datasets.")
    parser.add_argument("--input_file", type=str, required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save visualization results")
    parser.add_argument("--cluster_key", type=str, default="stlearn_louvain", help="Key for clustering in AnnData.obs (default: stlearn_louvain)")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    adata = sc.read_h5ad(args.input_file)

    base_name = os.path.basename(args.input_file)
    if "dlpfc" in base_name.lower():
        name = "dlpfc"
    elif "atac" in base_name.lower():
        name = "atac"
    else:
        name = "sp"
    # Main clustering (louvain)
    if args.cluster_key in adata.obs:
        y = adata.obs[args.cluster_key].astype(str)
        sizes, props = cluster_stats(y)
        plot_cluster_sizes(sizes, args.output_dir, f"{name}_{args.cluster_key}")
        plot_umap_from_obs(adata, args.cluster_key, args.output_dir, f"{name}_{args.cluster_key}")
        plot_spatial(adata, y, args.output_dir, f"{name}_{args.cluster_key}")
        # markers = marker_genes(adata, y, top_n=10)
        # plot_marker_heatmap(adata, y, markers, args.output_dir, f"{name}_{args.cluster_key}")
    else:
        print(f"Clustering key '{args.cluster_key}' not found in AnnData.obs.")
    # Kmeans clustering if present
    display_kmeans_clustering(adata, args.output_dir, name)
if __name__ == "__main__":
    main()
