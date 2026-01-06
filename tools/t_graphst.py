import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

from GraphST import GraphST
from GraphST.utils import clustering


import scipy.io as sio
from scipy.sparse import csr_matrix

import argparse

import dataloader as dl

def train_model(adata, device):
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    return adata

def cluster(adata, n_clusters, radius):
   
    clustering(adata, n_clusters, radius=radius, method='mclust', refinement=True)

    # rename obs: domain to obs:graphst_cluster
    adata.obs['graphst_cluster'] = adata.obs['domain']

    return adata



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    parser.add_argument("--input", type=str, required=True, help="Name of the input AnnData file.")
    parser.add_argument("--n_clusters", type=int, default=12, help="Number of clusters for clustering.")
    parser.add_argument("--radius", type=float, default=50, help="Radius for spatial graph construction.")
    args = parser.parse_args()

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    # Where R is located
    conda_env = "/opt/conda/envs/t_graphst/"

    os.environ['R_HOME'] = os.path.join(conda_env, "lib/R")
    os.environ['PATH'] = os.path.join(conda_env, "bin") + ":" + os.environ['PATH']



    adata = dl.open_anndata("../files/" + args.input)

    adata = train_model(adata, device)
    adata = cluster(adata, n_clusters=args.n_clusters, radius=args.radius)

    output = "../files/graphst_" + args.input

    adata.write_h5ad(output)