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
index_list = None  # reference order from the first CSV

for f in sorted(csv_files):
    df = pd.read_csv(os.path.join(data_path_dir, f), header=0, index_col=0)

    # ensure index is barcode strings
    df.index = df.index.astype(str)

    method_name = f.rsplit(".csv", 1)[0]
    df.columns = [method_name]

    # set reference barcode order from first file
    if index_list is None:
        index_list = df.index.tolist()
    else:
        # enforce same order + check for missing/extra barcodes
        if set(df.index) != set(index_list):
            missing = set(index_list) - set(df.index)
            extra = set(df.index) - set(index_list)
            raise ValueError(
                f"{f}: barcode set mismatch. missing={list(missing)[:10]}, extra={list(extra)[:10]}"
            )
        df = df.reindex(index_list)

    tools_df.append(df)

combined_df = pd.concat(tools_df, axis=1)

if index_list is None:
    raise ValueError(f"No CSV files found in {data_path_dir}")

combined_df = pd.concat(tools_df, axis=1)

# strict check: ensure reindex did not introduce NaNs
if combined_df.isna().any().any():
    bad_cols = combined_df.columns[combined_df.isna().any()].tolist()
    raise ValueError(f"NaNs after reindexing in columns: {bad_cols}")

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
#clsSim = nmi_similarity_rows(cluster_segments)
#print('finished computing cluster similarities: NMI\n')
clsSimRW = compute_pts_ii(clsSim, t=10)
print('finished computing cluster similarities: Random walk\n')
method_spot_mats = build_method_spot_matrices(bcs, clsSimRW)
print("Number of method spot matrices: ", len(method_spot_mats))
for idx, mat in enumerate(method_spot_mats):
    print(f"Method {methods[idx]} spot matrix shape: ", mat.shape)
H, w = ensdd_weighted_consensus(method_spot_mats, niter=100, epsilon=1e-5)
print("Consensus matrix H shape: ", H.shape)
print('finished computing consensus matrix\n')
print("Method weights: ", w)
H_csr = dense_similarity_to_csr_knn(H, k=50, tau=0.1, keep_diagonal=False, make_symmetric=True, dtype=np.float32)
print('finished building consensus kNN graph\n')
#labels, info = leiden_final(H_csr, K=K, seed=0, partition_type="rb")
#print('finished Leiden clustering with \n')
#print("Leiden info: ", info)

labels, info = leiden_overcluster_then_merge(A_csr=H_csr, H=H, K=K, seed=0, partition_type="rb", r_start=0.1,)
print('finished Leiden clustering with \n')
#print("Leiden info: ", info)

#labels = hc_from_consensus(H, K=K, linkage_method="average")
#print('finished Hierarchical Clustering to get final labels\n')


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