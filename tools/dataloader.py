import scanpy as sc
import argparse
import scipy.sparse as sp
import pandas as pd
from pathlib import Path



def load_data(file_path, file_type='h5ad'):
    """
    Load single-cell data from a specified file.

    Parameters:
    - file_path (str): The path to the data file.
    - file_type (str): The type of the data file. Default is 'h5ad'.

    Returns:
    - AnnData: The loaded single-cell data.
    """
    if file_type == 'h5ad':
        adata = sc.read_h5ad(file_path)
    elif file_type == 'loom':
        adata = sc.read_loom(file_path)
    elif file_type == 'mtx':
        adata = sc.read_mtx(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    return adata


def open_anndata(file_path):
    """
    Open an AnnData file.

    Parameters:
    - file_path (str): The path to the AnnData file.

    Returns:
    - AnnData: The opened AnnData object.
    """
    adata = sc.read_h5ad(file_path)
    return adata

def shared_preprocessing(adata):
    # Basic filtering
    sc.pp.filter_genes(adata, min_cells=3)

    # Library size normalization
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Highly variable genes (important!)
    sc.pp.highly_variable_genes(
        adata,
        flavor="cell_ranger",
        n_top_genes=3000
    )
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Scaling
    sc.pp.scale(adata, max_value=10)

    return adata

def save_cluster(output_path, adata, index_name):
    
    # check if directory exists for clusters, if not create it
    # Path("../files/" + user + "/clusters/").mkdir(parents=True, exist_ok=True)
    
    clusters = pd.DataFrame({
        "barcode": adata.obs_names,
        "spatial_cluster": adata.obs[index_name]
    })
    clusters.to_csv(output_path, index=False)
