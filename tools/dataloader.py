import scanpy as sc
import argparse


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

