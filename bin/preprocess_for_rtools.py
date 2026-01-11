import anndata
import pandas as pd
import numpy as np
import argparse as ap
import os


parser = ap.ArgumentParser(description="Preprocess AnnData for R tools")
parser.add_argument("--input", type=str, required=True, help="Input AnnData .h5ad file")
parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
args = parser.parse_args()

# If output directory does not exist, create it
os.makedirs(args.output_dir, exist_ok=True)

# Load your AnnData object
adata = anndata.read_h5ad(args.input)  # <-- change to your file


# Export counts (genes as rows, cells as columns)
counts = adata.X
if hasattr(counts, "toarray"):
    counts = counts.toarray()  # type: ignore
# adata.X is (cells, genes), so transpose to (genes, cells)
counts_df = pd.DataFrame(
    counts.T,  # type: ignore
    index=adata.var_names,
    columns=adata.obs_names
)

# Save in compressed format in output directory
counts_df_path = os.path.join(args.output_dir, "counts.csv.gz")
counts_df.to_csv(counts_df_path, compression="gzip")

# Export cell metadata (colData)
colData = adata.obs.copy()

colData["array_row"] = colData["centroid_y"]
colData["array_col"] = colData["centroid_x"]

colData_path = os.path.join(args.output_dir, "colData.csv")
colData.to_csv(colData_path)

# Export gene metadata (rowData)
rowData = adata.var.copy()
rowData_path = os.path.join(args.output_dir, "rowData.csv")
rowData.to_csv(rowData_path)