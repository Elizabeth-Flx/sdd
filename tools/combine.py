import dataloader as dl
import argparse
import scanpy as sc
import os
import pandas as pd

def add_cluster(adata, cluster_path, cluster_name):
    cluster = pd.read_csv(cluster_path)
    cluster = cluster.set_index("barcode")
    adata.obs[cluster_name] = cluster.loc[adata.obs_names, "spatial_cluster"]
    return adata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add clustering results to AnnData object.")
    parser.add_argument("--user", type=str, required=True, help="User name for file paths.")
    parser.add_argument('--clusters_names', type=list, nargs='*', default=[], help='List of cluster names to combine')
    args = parser.parse_args()

    adata = dl.open_anndata("../files/" + args.user + "/anndata.h5ad")
    
    # adata = dl.shared_preprocessing(adata)
    
    if len(args.clusters_names) == 0:
        # if no clusters provided use all in directory
        cluster_files = os.listdir("../files/" + args.user + "/clusters/")
        
        for file in cluster_files:
            print("Combining " + file)
            cluster_name = file.split(".csv")[0]
            adata = add_cluster(adata, "../files/" + args.user + "/clusters/" + file, cluster_name)
    else:
        for cluster_name in args.clusters_names:
            print("Combining " + cluster_name)
            cluster_path = args.user + "/clusters/" + cluster_name + ".csv"
            adata = add_cluster(adata, "../files/" + cluster_path, cluster_name)
        
    adata.write_h5ad("../files/" + args.user + "/anndata_combined.h5ad")
