#! /usr/bin/env python3

import scanpy as sc
import os
import pandas as pd
import argparse as ap

# # List of input files
# input_files = [
#     "data/tangram_output/sp_rna.h5ad",
#     "data/tangram_output/sp_atac.h5ad"
# ]

def convert_x_to_csv(filename: str, output_dir: str):
    adata = sc.read_h5ad(filename)
    # Convert X to dense if needed
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray() # type: ignore
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    outname = os.path.splitext(filename)[0] + "_counts.csv"
    outname = os.path.join(output_dir, os.path.basename(outname))
    df.to_csv(outname)
    # print(f"Written: {outname}")



if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Convert AnnData X matrices to CSV")
    parser.add_argument("--input_files", required=True, nargs='+', help="List of input AnnData .h5ad files")
    parser.add_argument("--output_dir", required=False, help="Output directory for CSV files", default=".")

    args = parser.parse_args()
    input_files = args.input_files
    output_dir = args.output_dir
    for infile in input_files:
        convert_x_to_csv(infile, output_dir)