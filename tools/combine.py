from . import dataloader as dl
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
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--csv_dir", type=str, required=True, help="Input csv directory")
    parser.add_argument("--output", type=str, required=True, help="Output anndata file")   
    
    parser.add_argument("--clusters_names", nargs='*', default=[], help="List of cluster names to add (without .csv). If empty, all clusters in the directory will be added.")
     
    args = parser.parse_args()

    adata = dl.open_anndata(args.input)
    
    # adata = dl.shared_preprocessing(adata)
    
    if len(args.clusters_names) == 0:
        # if no clusters provided use all in directory
        cluster_files = os.listdir(args.csv_dir)
        
        for file in cluster_files:
            if file.endswith(".csv") == False:
                continue
            print("Combining " + file)
            cluster_name = file.split(".csv")[0]
            adata = add_cluster(adata, args.csv_dir + "/" + file, cluster_name)
    else:
        for cluster_name in args.clusters_names:
            print("Combining " + cluster_name)
            cluster_path = args.csv_dir + "/" + cluster_name + ".csv"
            adata = add_cluster(adata, cluster_path, cluster_name)

    adata.write_h5ad(args.output)