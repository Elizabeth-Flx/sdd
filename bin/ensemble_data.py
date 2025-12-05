#! /usr/bin/env python3

# NOTE: This is still a dummy script, changes will likely be necessary (exact formats etc. are still unknown)

import os
import pandas as pd
import numpy as np
import argparse as ap
import scanpy as sc


def create_matrix_for_ensemble(sp_list, tool_list, annotation_key="celltype"):
    """
    Create a matrix for ensemble methods from the given list of spatial data.
    @param sp_list: List of spatial data objects
    @return: Matrix cells*tools
    """
    # 2dimensional matrix
    
   # Get the set of all unique cell names (union across all tools)
    all_cells = set()
    for ad in sp_list:
        all_cells.update(ad.obs_names)
    all_cells = sorted(all_cells)

    # Build a DataFrame: rows=cells, columns=tools
    data = {}
    for i, ad in enumerate(sp_list):
        # Map cell name to annotation, fill missing with np.nan
        ann = ad.obs[annotation_key] if annotation_key in ad.obs else pd.Series(np.nan, index=ad.obs_names)
        data[tool_list[i]] = pd.Series(ann.values, index=ad.obs_names)
    df = pd.DataFrame(data, index=all_cells)
    return df


if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Create a matrix for ensemble methods from spatial data.")
    parser.add_argument('--input_dir', help='Directory containing spatial data objects from tangram')
    parser.add_argument('--output_file', help='Output file to save the matrix')
    args = parser.parse_args()

    # Iterate through directory, each file has the format tangram_$tool.h5ad
    sp_list = []
    tool_list = []

    for sp_data_file in os.listdir(args.input_dir):
        # file path could start with path/to/file
        if "tangram_" in sp_data_file and sp_data_file.endswith(".h5ad"):
            tool_name = sp_data_file.split("tangram_")[-1].replace(".h5ad", "")
            sp_list.append(sc.read_h5ad(os.path.join(args.input_dir, sp_data_file)))
            tool_list.append(tool_name)

    matrix = create_matrix_for_ensemble(sp_list, tool_list)

    # Write matrix to file
    matrix.to_csv(args.output_file)