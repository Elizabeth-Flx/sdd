import os
import torch
from pathlib import Path
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

from GraphST import GraphST
from GraphST.utils import clustering


import scipy.io as sio
from scipy.sparse import csr_matrix

import argparse

from . import dataloader as dl

def train_model(adata, device):
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    return adata

def cluster(adata, 
            cluster_method, 
            k_clusters, 
            resolution, 
            radius):
    
    if cluster_method == 'mclust':
        clustering(adata, k_clusters, radius=radius, method='mclust', refinement=True)
    elif cluster_method == 'leiden':
        sc.pp.neighbors(adata, use_rep="emb")
        sc.tl.leiden(adata, resolution=resolution, key_added="domain")
    elif cluster_method == 'louvain':
        sc.pp.neighbors(adata, use_rep="emb")
        sc.tl.louvain(adata, resolution=resolution, key_added="domain")

    
    
    
    # rename obs: domain to obs:graphst_cluster
    adata.obs['graphst_cluster'] = adata.obs['domain'].astype("category")

    return adata



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--output", type=str, required=True, help="Output csv file")    
    
    # parser.add_argument("--input", type=str, required=True, help="Name of the input AnnData file.")
    parser.add_argument("--cluster_method", type=str, default="mclust", help="Clustering method: mclust, leiden, louvain.")
    parser.add_argument("--k_clusters", type=int, default=12, help="Number of clusters for clustering (used in mclust).")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for clustering (used in leiden and louvain).")
    parser.add_argument("--radius", type=float, default=50, help="Radius for spatial graph construction.")
    parser.add_argument("--do_normlog", type=bool, default=True, help="Whether to perform preprocessing.")

    args = parser.parse_args()

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    # Where R is located
    conda_env = "/opt/conda/envs/t_graphst/"

    os.environ['R_HOME'] = os.path.join(conda_env, "lib/R")
    os.environ['PATH'] = os.path.join(conda_env, "bin") + ":" + os.environ['PATH']

    adata = dl.open_anndata(args.input)
    
    adata = dl.shared_preprocessing(adata, args.do_normlog)

    adata = train_model(adata, device)
    adata = cluster(adata, cluster_method=args.cluster_method, k_clusters=args.k_clusters, resolution=args.resolution, radius=args.radius)

    dl.save_cluster(args.output, adata, "graphst_cluster")
