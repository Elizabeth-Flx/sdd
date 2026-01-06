#! /usr/bin/env python3

import scanpy as sc
import argparse
from GraphST import GraphST
from GraphST.utils import clustering


def main(args):
    # Load AnnData object from h5ad file
    adata = sc.read_h5ad(args.input_h5ad)

    # Set up device
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    device = "cpu"
    
    # Run GraphST
    model = GraphST.GraphST(
        adata,
        device,
    )
    adata = model.train()

    # Clustering
    clustering(
        adata,
        args.n_clusters,
        radius=args.radius,
        method=args.method,
        refinement=True
    )

    # Optionally, save the result
    if args.output_h5ad:
        adata.write_h5ad(args.output_h5ad)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GraphST on an existing AnnData object (h5ad file).")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to input AnnData h5ad file")
    parser.add_argument("--output_h5ad", type=str, required=False, help="Path to save output AnnData h5ad file")
    parser.add_argument("--n_clusters", type=int, default=10, help="Number of clusters for clustering")
    parser.add_argument("--radius", type=float, default=150, help="Radius for clustering (in micrometers)")
    parser.add_argument("--method", type=str, default="mclust", help="Clustering method (default: mclust)")
    args = parser.parse_args()
    main(args)
