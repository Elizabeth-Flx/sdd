import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse as ap
import seaborn as sns
# Add cluster colorbar
from matplotlib.colors import ListedColormap


def get_cluster_labels(adata):
    # Prefer mclust if available
    if "mclust" in adata.obs:
        return adata.obs["mclust"].astype(str)
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


def plot_umap_from_obs(adata, out_path, name):
    out_file = f"{out_path}/umap_{name}.png"
    # print(f"Saving UMAP plot to {out_file}")
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)  # pyright: ignore[reportCallIssue]
    figs = sc.pl.umap(adata, color='mclust', show=False, return_fig=True)
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
        # print(f"Cluster {cl} top markers: {len(top_genes)} genes: {top_genes}")
        markers[cl] = set(top_genes)
    return markers


def plot_marker_heatmap(adata, y, markers, out_path, name):
    all_markers = list(set.union(*markers.values()))
    sc.pl.heatmap(adata, all_markers, y.name) # pyright: ignore[reportCallIssue]
    plt.savefig(f"{out_path}/marker_heatmap_{name}.png")
    plt.close()



def plot_cell_gene_heatmap(
    adata, 
    cluster_labels, 
    genes=None, 
    out_path=".", 
    name="heatmap", 
    cluster_rows=True, 
    cluster_cols=True, 
    max_cells=1000,
    max_genes=50
):
    """
    Plots a heatmap where rows are cells (with cluster annotation) and columns are genes.
    Optionally clusters rows and/or columns.
    - adata: AnnData object
    - cluster_labels: pd.Series or list of cluster labels (length = n_cells)
    - genes: list of genes to plot (default: all genes)
    - out_path: directory to save the plot
    - name: filename suffix
    - cluster_rows: whether to cluster rows (cells)
    - cluster_cols: whether to cluster columns (genes)
    - max_cells: maximum number of cells to plot (for performance)
    """


    # Subset genes if provided, else use top highly variable genes (default: max_genes)
    if genes is None:
        # Try to use highly variable genes if present
        if 'highly_variable' in adata.var.columns:
            genes = adata.var_names[adata.var['highly_variable']].tolist()
            if len(genes) > max_genes:
                genes = genes[:max_genes]
        else:
            # If not present, compute variance and select top max_genes
            X = adata.X
            if hasattr(X, "toarray"):
                X = X.toarray()
            gene_vars = np.var(X, axis=0)
            top_idx = np.argsort(gene_vars)[-max_genes:]
            genes = adata.var_names[top_idx].tolist()
    # Subset cells if too many
    if adata.n_obs > max_cells:
        idx = np.random.choice(adata.n_obs, max_cells, replace=False)
    else:
        idx = np.arange(adata.n_obs)
    # Extract expression matrix
    X = adata[idx, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    # Build DataFrame
    df = pd.DataFrame(X, columns=genes)
    df['cluster'] = pd.Series(cluster_labels).iloc[idx].values
    # Sort by cluster for better visualization
    df = df.sort_values('cluster')
    clusters = df['cluster']
    df = df.drop(columns='cluster')
    # Create color mapping for clusters
    unique_clusters = pd.Series(clusters).unique()
    lut = dict(zip(unique_clusters, sns.color_palette("tab20", len(unique_clusters))))
    row_colors = clusters.map(lut)
    # Plot with row_colors
    g = sns.clustermap(
        df,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        row_colors=row_colors,
        yticklabels=False,
        xticklabels=True,
        cmap="viridis",
        figsize=(12, 8),
        dendrogram_ratio=(.1, .2),
        cbar_pos=(0.05, 0.2, 0.03, 0.4) 
    )
    # Add legend for clusters
    for label in unique_clusters:
        g.ax_row_colors.bar(0, 0, color=lut[label], label=str(label), linewidth=0)
    g.ax_row_colors.legend(
        title="Cluster",
        bbox_to_anchor=(1.15, 0.5),
        loc='center left',
        ncol=1,
        frameon=True,
        borderaxespad=0.0
    )
    g.ax_heatmap.set_ylabel("Cells (clustered)")
    g.ax_heatmap.set_xlabel("Genes")
    # After creating the clustermap and before saving:
    plt.title("Expr.", y=1.0, loc = "center")
    # g.fig.subplots_adjust(top=0.92)  # Adjust as needed for your layout
    plt.savefig(os.path.join(out_path, f"cell_gene_heatmap_{name}.png"))
    plt.close()


def main():
    parser = ap.ArgumentParser(description="Compare clustering results between two AnnData datasets.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing input AnnData files")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save comparison results")

    args = parser.parse_args()

    # Make output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Use tool name in filenames
    dlpfc_file = os.path.join(args.input_dir, "sp_dlpfc_graphst.h5ad")
    atac_file = os.path.join(args.input_dir, "sp_atac_graphst.h5ad")
    sp_file = os.path.join(args.input_dir, "sp_graphst.h5ad")

    ad_dlpfc = sc.read_h5ad(dlpfc_file)
    ad_atac = sc.read_h5ad(atac_file)
    ad_sp = sc.read_h5ad(sp_file)

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

    # Cell-gene heatmaps (all genes, default clustering)
    plot_cell_gene_heatmap(ad_dlpfc, y_dlpfc, out_path=args.output_dir, name="dlpfc", cluster_rows=False, cluster_cols=True, max_genes=50)
    plot_cell_gene_heatmap(ad_atac, y_atac, out_path=args.output_dir, name="atac", cluster_rows=False, cluster_cols=True, max_genes=50)
    plot_cell_gene_heatmap(ad_sp, y_sp, out_path=args.output_dir, name="sp", cluster_rows=False, cluster_cols=True, max_genes=50)


if __name__ == "__main__":
    main()