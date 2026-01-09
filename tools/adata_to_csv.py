import scanpy as sc
import pandas as pd
import scipy.sparse as sp


def create_csvs(adata_path, csv_dir="../files/csv/"):
    adata = sc.read_h5ad(adata_path)
    X = adata.X
    
    if sp.issparse(X):
        countData = pd.DataFrame.sparse.from_spmatrix(
            X.T,
            index=adata.var_names,
            columns=adata.obs_names
        )
    else:
        countData = pd.DataFrame(
            X.T,
            index=adata.var_names,
            columns=adata.obs_names
        )
        
    rowData = adata.var.copy()
    rowData.insert(0, "gene_id", rowData.index)

    colData = adata.obs.copy()

    countData.to_csv(csv_dir + "/counts.csv.gz", compression="gzip")
    rowData.to_csv(csv_dir + "/rowData.csv")
    colData.to_csv(csv_dir + "/colData.csv")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert AnnData to CSV files for R-based tools.")
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for CSV files")    
    
    args = parser.parse_args()

    create_csvs(args.input, args.output_dir)