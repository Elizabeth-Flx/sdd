import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

import networkx as nx
import community as community_louvain
from ensemble_helper import *
import os
import argparse

parser = argparse.ArgumentParser(description="Generating ensemble domain detection clusters.")

parser.add_argument("--input_csvs", type=str, required=True, help="Input tool csv files directory")
parser.add_argument('--cluster_k', type=int, required=False, default=12, help='Number of clusters for final clustering (default: 12)')
args = parser.parse_args()
data_path_dir = args.input_csvs
K = args.cluster_k
print('Input csv files directory:', data_path_dir)

# for all csvs in dir: read in each csv (only one column per csv) to dataframe to use
import os
csv_files = [f for f in os.listdir(data_path_dir) if f.endswith('.csv')]
tools_df = []
# see if len are identical
indices = []
lengths = []
for f in csv_files:
    df = pd.read_csv(os.path.join(data_path_dir, f), header=0, index_col=0)
    indices.append(df.index.tolist())
    lengths.append(df.shape[0])
    method_name = f.split('.csv')[0]
    df.columns = [method_name]
    # set index first to 1: so index starting from 1
    df.index = range(1, df.shape[0] + 1)
    tools_df.append(df)
    
combined_df = pd.concat(tools_df, axis=1)
if len(set(lengths)) != 1:
    raise ValueError("CSV files have different number of rows!")

# check if all indices identical
if not all(x == indices[0] for x in indices):
    raise ValueError("CSV files have different indices!")
index_list = indices[0]

print("Combined dataframe shape: ", combined_df.shape)
print("Combined dataframe tools: ", combined_df.columns.tolist())
print('Finished reading in csv files\n')

# use all predefined methods
print('Starting Ensemble calculation...')
baseCls, spot_ids, methods = encode_cluster_labels(combined_df)
print('methods: ', methods)
print('finished encoding cluster labels\n')
bcs, nClsOrig, C_total = make_global_cluster_ids(baseCls)
print('finished making global cluster ids\n')
cluster_segments = build_cluster_segments(bcs, C_total)
print('finished building cluster segments\n')
clsSim = jaccard_similarity_rows(cluster_segments)
print('finished computing cluster similarities: Jaccard\n')
clsSimRW = compute_pts_ii(clsSim, t=20)
print('finished computing cluster similarities: Random walk\n')
method_spot_mats = build_method_spot_matrices(bcs, clsSimRW)
print("Number of method spot matrices: ", len(method_spot_mats))
for idx, mat in enumerate(method_spot_mats):
    print(f"Method {methods[idx]} spot matrix shape: ", mat.shape)
H, w = ensdd_weighted_consensus(method_spot_mats, niter=100, epsilon=1e-5)
print("Consensus matrix H shape: ", H.shape)
print('finished computing consensus matrix\n')
print("Method weights: ", w)
H_csr = dense_similarity_to_csr_knn(H, k=5, keep_diagonal=False, make_symmetric=True, dtype=np.float32)
print('finished building consensus kNN graph\n')
labels, info = leiden_final(H_csr, K=K, seed=0, partition_type="rb")
print('finished Leiden clustering with \n')
#print("Leiden info: ", info)

# to output dataframe format
leiden_df = pd.DataFrame(
    {"spatial_cluster": labels},
    index=index_list
)
leiden_df.index.name = "barcode"

# save to csv in same dir as input csvs
path_ensemble = os.path.join(data_path_dir, 'ensemble.csv')
leiden_df.to_csv(path_ensemble)

print(f"Ensemble clustering labels saved to {path_ensemble}")
print('Ensemble clustering process finished.')