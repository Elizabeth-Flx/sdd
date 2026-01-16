import scanpy as sc
import scanit
import numpy as np
import argparse
from sklearn.cluster import KMeans

import torch
from . import dataloader as dl

def train_model(adata, device=torch.device('cpu')):
    scanit.tl.spatial_graph(
        adata,
        method="alpha shape",
        alpha_n_layer=2,
        knn_n_neighbors=15
    )
        
    if isinstance(adata.X, (np.ndarray)):
        X_dense = adata.X
    else:
        X_dense = adata.X.toarray()

    adata.X = X_dense

    scanit.tl.spatial_representation(
        adata,
        n_h=30,
        n_epoch=1000,
        lr=0.001,
        device=device,
        n_consensus=1,
        projection="mds",
        python_seed=0,
        torch_seed=0,
        numpy_seed=0
    )
    return adata

def cluster(adata, 
            cluster_method="leiden", 
            k_clusters=12,
            resolution=1.0,
            n_neighbors=15):

    if cluster_method == "kmeans":
        kmeans = KMeans(n_clusters=k_clusters, random_state=0)
        adata.obs["scanit_cluster"] = kmeans.fit_predict(adata.obsm["X_scanit"]).astype(str)
    elif cluster_method == "louvain":
        sc.pp.neighbors(adata, use_rep="X_scanit", n_neighbors=n_neighbors)
        sc.tl.louvain(adata, resolution=resolution, key_added="scanit_cluster")
    elif cluster_method == "leiden":
        sc.pp.neighbors(adata, use_rep="X_scanit", n_neighbors=n_neighbors)
        sc.tl.leiden(adata, resolution=resolution, key_added="scanit_cluster")

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--output", type=str, required=True, help="Output csv file")    
    
    parser.add_argument("--cluster_method", type=str, default="kmeans", help="Clustering method: kmeans, leiden, louvain.")
    parser.add_argument("--k_clusters", type=int, default=12, help="Number of clusters for kmeans.")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for Louvain clustering.")
    parser.add_argument("--n_neighbors", type=int, default=19, help="Number of neighbors for graph construction.")
    parser.add_argument("--do_normlog", type=bool, default=True, help="Whether to perform preprocessing.")

    args = parser.parse_args()

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    adata = dl.open_anndata(args.input)

    adata = dl.shared_preprocessing(adata, args.do_normlog)
    
    adata = train_model(adata, device)
    adata = cluster(adata, cluster_method=args.cluster_method, k_clusters=args.k_clusters, resolution=args.resolution, n_neighbors=args.n_neighbors)

    dl.save_cluster(args.output, adata, "scanit_cluster")
