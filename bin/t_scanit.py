import scanpy as sc
import scanit
import numpy as np
import argparse

import torch

def preprocessing(adata):
    # Preprocessing
    sc.pp.normalize_total(adata)  # type: ignore
    sc.pp.log1p(adata)  # type: ignore
    sc.pp.scale(adata)  # type: ignore

    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)  # type: ignore
    adata = adata[:, adata.var["highly_variable"]]
    
    return adata

def train_model(adata, device=torch.device('cpu')):
    scanit.tl.spatial_graph(
        adata,
        method="alpha shape",
        alpha_n_layer=2,
        knn_n_neighbors=15
    )

    scanit.tl.spatial_representation(
        adata,
        n_h=30,
        n_epoch=2000,
        lr=0.001,
        device=device,
        n_consensus=1,
        projection="mds",
        python_seed=0,
        torch_seed=0,
        numpy_seed=0
    )
    return adata

def cluster(adata,  louvain_resolution=1.19):

    # Clustering with louvain
    sc.pp.neighbors(adata, use_rep="X_scanit", n_neighbors=15)
    sc.tl.leiden(adata, resolution=louvain_resolution, key_added="scanit_louvain")

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to input AnnData h5ad file")
    parser.add_argument("--output_h5ad", type=str, required=False, help="Path to save output AnnData h5ad file")
    parser.add_argument("--louvain_resolution", type=float, default=1.19, help="Resolution parameter for Louvain clustering.")
    args = parser.parse_args()

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')


    adata = sc.read_h5ad(args.input_h5ad)

    adata = preprocessing(adata)
    adata = train_model(adata, device)
    adata = cluster(adata, louvain_resolution=args.louvain_resolution)

    output = args.output_h5ad

    adata.write_h5ad(output)