#! /usr/bin/env python3
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse as ap
import seaborn as sns
from matplotlib.colors import ListedColormap  # type: ignore


def get_cluster_labels(adata, cluster_key):
    if cluster_key in adata.obs:
        return adata.obs[cluster_key].astype(str)
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
    plt.scatter(x, y_coord, c=cluster_codes, cmap='tab20', s=7, alpha=0.8)
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


def plot_umap_from_obs(adata, out_path, name, cluster_key):
    out_file: str = f"{out_path}/umap_{name}.png"
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)  # pyright: ignore[reportCallIssue]
    figs = sc.pl.umap(adata, color=str(cluster_key), show=False, return_fig=True)
    if isinstance(figs, list):
        fig = figs[0]
    else:
        fig = figs
    fig.savefig(out_file)  # type: ignore
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
        g.ax_row_colors.bar(0, 0, color=lut[label], label=str(label), linewidth=0)  # type: ignore
        
    g.ax_row_colors.legend(  # type: ignore
        title="Cluster",
        bbox_to_anchor=(1.15, 0.5),
        loc='center left',
        ncol=1,
        frameon=True,
        borderaxespad=0.0
    )  
    g.ax_heatmap.set_ylabel("Cells (clustered)")
    g.ax_heatmap.set_xlabel("Genes")
    plt.title("Expr.", y=1.0, loc = "center")
    # g.fig.subplots_adjust(top=0.92)  # Adjust as needed
    plt.savefig(os.path.join(out_path, f"cell_gene_heatmap_{name}.png"))
    plt.close()




def main():
    parser = ap.ArgumentParser(description="Visualize clusters from CSV and h5ad file.")
    parser.add_argument('--csv', required=True, help='Path to cluster CSV file (first column: obs_names, second: cluster label)')
    parser.add_argument('--h5ad', required=True, help='Path to original h5ad file')
    parser.add_argument("--tool_name", type=str, required=True, help="Name of the tool (used in filenames)")
    parser.add_argument('--output_dir', required=True, help='Output directory for figures')
    args = parser.parse_args()

    # Make output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    if "atac" in args.h5ad.lower():
        name = "atac"
    else:
        name = "sp"

    # Load cluster assignments from CSV
    cluster_df = pd.read_csv(args.csv)
    if cluster_df.shape[1] < 2:
        raise ValueError("CSV must have at least two columns: obs_names and cluster label.")
    obs_names = cluster_df.iloc[:,0].astype(str)
    clusters = cluster_df.iloc[:,1].astype(str)
    cluster_map = dict(zip(obs_names, clusters))

    # Load AnnData
    adata = anndata.read_h5ad(args.h5ad)

    # Match and add cluster info to AnnData.obs
    adata.obs['cluster'] = adata.obs_names.map(lambda x: cluster_map.get(str(x), None))  # type: ignore
    missing = adata.obs['cluster'].isnull().sum()
    if missing > 0:
        print(f"Warning: {missing} cells in h5ad not found in CSV. They will have cluster=None.")

    y_sp = adata.obs['cluster']

    # UMAP plots
    plot_umap_from_obs(adata, args.output_dir, name, "cluster")

    # Spatial plots
    plot_spatial(adata, y_sp, args.output_dir, name)

    # Cell-gene heatmaps (top 50 genes)
    plot_cell_gene_heatmap(adata, y_sp, out_path=args.output_dir, name=name, cluster_rows=False, cluster_cols=True, max_genes=50)


if __name__ == '__main__':
	main()
