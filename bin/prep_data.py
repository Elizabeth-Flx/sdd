#! /usr/bin/env python3

import scanpy as sc
import pandas as pd 
from scipy import io
import os
import anndata as ad

DATAPATH = "/nfs/data3/mopitas/sddproj/data/dlpfc"

X = io.mmread(os.path.join(DATAPATH, "DLPFC_151671_rawcounts.mtx"))       # reads as COO sparse matrix
X = X.tocsr() 
obs_data = pd.read_csv(os.path.join(DATAPATH, "DLPFC_151671_cell_data.csv"), index_col=0)
var_data = pd.read_csv(os.path.join(DATAPATH, "DLPFC_151671_gene_data.csv"), index_col=0)
spatial_coords = pd.read_csv(os.path.join(DATAPATH, "DLPFC_151671_spatial_coords.csv"), index_col=0)

adata = ad.AnnData(X=X.T, obs=obs_data, var=var_data)
adata.obsm['spatial'] = spatial_coords.values
adata.var_names = adata.var['gene_name']
adata.var_names_make_unique()
if 'gene_name' in adata.var.columns:
    adata.var = adata.var.drop(columns=['gene_name'])
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=500)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.pca(adata)

output_path = "data/"
# Save processed data as h5ad file
adata.write_h5ad(os.path.join(output_path, "DLPFC_151671_processed.h5ad"))