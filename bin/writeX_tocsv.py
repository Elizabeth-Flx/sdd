import scanpy as sc
import os
import pandas as pd

# List of input files
input_files = [
    "data/tangram_output/sp_rna.h5ad",
    "data/tangram_output/sp_atac.h5ad"
]



for infile in input_files:
    adata = sc.read_h5ad(infile)
    # Convert X to dense if needed
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray() # type: ignore
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    outname = os.path.splitext(infile)[0] + "_counts.csv"
    df.to_csv(outname)
    print(f"Written: {outname}")
