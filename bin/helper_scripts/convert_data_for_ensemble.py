#! /usr/bin/env python3

# Import necessary libraries
import scanpy as sc
import pandas as pd

# Load integrated data for each tool (some csv, some h5ad) and extract cluster labels
# Write a csv file for each tool in new directory, with a sinle column containing cluster labels

# The column is named spatial_cluster in the csv files -> bass & bayesspace
# For the other tools, the cluster keys are as follows, stored in a .h5ad file:
# TOOL_CLUSTER_KEYS["graphst"]="mclust"
# TOOL_CLUSTER_KEYS["scanit"]="scanit_louvain"
# TOOL_CLUSTER_KEYS["stlearn"]="stlearn_louvain"


TOOL_CLUSTER_KEYS = {
    "stlearn": "stlearn_louvain",
    # "graphst": "mclust",
    # "scanit": "scanit_louvain",
    # "bass": "spatial_cluster",
    # "bayesspace": "spatial_cluster"
}

def get_cluster_label_csv (csv_path, cluster_key, out_path):
    df = pd.read_csv(csv_path, index_col=0)
    if cluster_key not in df.columns:
        raise ValueError(f"Cluster key {cluster_key} not found in {csv_path}.")
    cluster_labels = df[cluster_key]
    out_df = pd.DataFrame({ "spatial_cluster": cluster_labels })
    # Write without column name or index or header
    out_df.to_csv(out_path, index=False, header=False)

def get_cluster_label_h5ad (h5ad_path, cluster_key, out_path):
    ad = sc.read_h5ad(h5ad_path)
    if cluster_key not in ad.obs:
        raise ValueError(f"Cluster key {cluster_key} not found in {h5ad_path}.")
    cluster_labels = ad.obs[cluster_key]
    out_df = pd.DataFrame({ "spatial_cluster": cluster_labels })
    # Write without column name or index or header
    out_df.to_csv(out_path, index=False, header=False)


def convert_data_for_ensemble(input_dir, output_dir):
    # Each tool has its own directory, directory name = tool name
    # There we have a file starting with sp_integrated_{tool}.csv or sp_integrated_rna.h5ad
    import os
    os.makedirs(output_dir, exist_ok=True)
    for tool, cluster_key in TOOL_CLUSTER_KEYS.items():
        tool_input_dir = os.path.join(input_dir, tool)
        if not os.path.isdir(tool_input_dir):
            print(f"Warning: Tool input directory {tool_input_dir} does not exist. Skipping.")
            continue
        # Find the relevant file
        csv_file = os.path.join(tool_input_dir, f"sp_integrated_{tool}.csv")
        h5ad_file = os.path.join(tool_input_dir, f"sp_integrated_rna.h5ad")
        out_file = os.path.join(output_dir, f"{tool}.csv")
        if os.path.isfile(csv_file):
            get_cluster_label_csv(csv_file, cluster_key, out_file)
            print(f"Wrote cluster labels for {tool} from CSV to {out_file}.")
        elif os.path.isfile(h5ad_file):
            get_cluster_label_h5ad(h5ad_file, cluster_key, out_file)
            print(f"Wrote cluster labels for {tool} from h5ad to {out_file}.")
        else:
            print(f"Warning: No input file found for tool {tool} in {tool_input_dir}. Skipping.")


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="Convert tool outputs to cluster label CSVs for ensemble.")
    parser.add_argument('--input_dir', required=True, help='Directory with tool output subdirectories')
    parser.add_argument('--output_dir', required=True, help='Directory to save cluster label CSVs')
    args = parser.parse_args()

    convert_data_for_ensemble(args.input_dir, args.output_dir)