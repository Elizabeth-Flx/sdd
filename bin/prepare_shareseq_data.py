#! /usr/bin/env python3

import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
import muon as mu
from muon import atac
from pathlib import Path
import sys
import re
import argparse


def get_rna_anndata(df_rna_counts, df_celltype, log_file=None):
    df_rna_counts.columns = df_rna_counts.columns.str.replace(',', '.', regex=False)
    df_celltype = df_celltype.set_index('rna.bc')
    common_cell_ids = df_celltype.index.tolist()
    df_rna_subset = df_rna_counts.reindex(columns=common_cell_ids).dropna(axis=1)
    # Verify that the columns now exactly match the common cell IDs
    if list(df_rna_subset.columns) != common_cell_ids:
        if log_file:
            with open(log_file, "a") as f:
                f.write("Error: Subsetting RNA counts failed to exactly match common cell IDs!\n")
        sys.exit(1)
        # Handle error - critical for alignment
    df_rna_transposed = df_rna_subset.T.copy()

    adata_rna = ad.AnnData(
        X=df_rna_transposed,
        obs=pd.DataFrame(index=df_rna_transposed.index), # Use transposed index (common cell IDs) for obs names
        var=pd.DataFrame(index=df_rna_transposed.columns) # Use transposed columns (gene IDs) for var names
    )
    adata_rna.obs = adata_rna.obs.merge(df_celltype[['celltype']], left_index=True, right_index=True, how='left')
    
    return adata_rna


def filter_rna_anndata(adata_rna, min_cells=2, min_genes=100, n_top_genes = 2000, flavor="seurat_v3", qc_metrics=True, normalize=True, log=True, highly_variable=True, pca = True, log_file=None):

    mito_count = adata_rna.var_names.str.startswith("Mt-").sum()
    # Write mito count to log file
    if log_file:
        with open(log_file, "a") as f:
            f.write("FILTERING RNA ANNDATA\n")
            f.write(f"Mitochondrial genes count: {mito_count}\n")
            # Write parameters used for filtering
            f.write(f"min_cells: {min_cells}\n\
                        min_genes: {min_genes}\n\
                        calculate qc_metrics: {qc_metrics}\n\
                        normalize counts: {normalize}\n\
                        logarithm (log10): {log}\n\
                        highly_variable: {highly_variable}\n")

    if qc_metrics:
        sc.pp.calculate_qc_metrics(
        adata_rna, inplace=True, log1p=True
    )

    # Filter data
    sc.pp.filter_cells(adata_rna, min_genes=min_genes) # pyright: ignore[reportCallIssue]
    sc.pp.filter_genes(adata_rna, min_cells=min_cells) # pyright: ignore[reportCallIssue]

    # Saving count data
    adata_rna.layers["counts"] = adata_rna.X.copy()

    if normalize:
        sc.pp.normalize_total(adata_rna) # pyright: ignore[reportCallIssue]
    if log:
        sc.pp.log1p(adata_rna)
    if highly_variable:
        sc.pp.highly_variable_genes(adata_rna, n_top_genes=n_top_genes, flavor=flavor) # pyright: ignore[reportCallIssue]
    if pca:
        sc.tl.pca(adata_rna)

    return adata_rna



def process_atac_data(celltype, atac, barcodes, peaks, ):
    atac.obs["barcode"] = barcodes[0].values
    peaks["peak_id"] = peaks["chrom"] + ":" + peaks["start"].astype(str) + "-" + peaks["end"].astype(str)
    peaks.set_index("peak_id", inplace=True)
    atac.var = peaks
    # Map RNA and ATAC barcodes
    rna_to_atac= dict(zip(celltype["rna.bc"], celltype["atac.bc"]))

    # Add ATAC barcode column based on RNA barcodes
    atac.obs["atac_barcode"] = [rna_to_atac.get(bc, "unknown") for bc in atac.obs["barcode"]]

    # Format ATAC barcodes for fragments file
    formatted_barcodes = []
    for bc in atac.obs["atac_barcode"]:
        components = re.findall(r'(R\d+\.\d+|P\d+\.\d+)', bc)
        formatted_bc = ','.join(components)
        formatted_barcodes.append(formatted_bc)
    
    # Add formatted barcodes
    atac.obs["formatted_barcode"] = formatted_barcodes

    # Set formatted barcodes as observation names
    atac.obs_names = atac.obs["formatted_barcode"].values

    # Add cell type information directly
    atac.obs["celltype"] = atac.obs["barcode"].map(dict(zip(celltype["rna.bc"], celltype["celltype"])))

    return atac


def add_metrics_atac(atac, fragment_min = 165, fragment_max = 12000, nucleosome_signal_threshold=2, nucleosome_signal = 2000000,  log_file=None):
    # Calculate general qc metrics using scanpy
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

    # Rename columns
    atac.obs.rename(
        columns={
            "n_genes_by_counts": "n_features_per_cell",
            "total_counts": "total_fragment_counts",
        },
        inplace=True,
    )
    # log-transform total counts and add as column
    atac.obs["log_total_fragment_counts"] = np.log10(atac.obs["total_fragment_counts"])
    # Calculate the nucleosome signal across cells
    # set n=10e3*atac.n_obs for rough estimate but faster run time
    mu.atac.tl.nucleosome_signal(atac, n=nucleosome_signal)

    # Add group labels for above and below the nucleosome signal threshold
    atac.obs["nuc_signal_filter"] = [
        "NS_FAIL" if ns > nucleosome_signal_threshold else "NS_PASS"
        for ns in atac.obs["nucleosome_signal"]
    ]

    if log_file is not None:
        with open(log_file, "a") as f:
            f.write("Nucleosome signal filter value counts:\n")
            f.write(str(atac.obs["nuc_signal_filter"].value_counts()))
            f.write("\n")
    pp_data = atac.copy()

    if log_file is not None:
        with open(log_file, "a") as f:
            f.write(f"Total number of cells: {pp_data.n_obs}\n")
            
    mu.pp.filter_obs(
        pp_data,
        "total_fragment_counts",
        lambda x:  (x >= fragment_min) &(x <= fragment_max),  # parameter for filtering
    )
    if log_file is not None:
        with open(log_file, "a") as f:
            f.write(f"Number of cells after filtering on total_fragment_counts: {pp_data.n_obs}\n")
            
    return pp_data





if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Prepare SHARE-seq data for downstream analysis.")
    parser.add_argument('--input_dir', required=True, help='Directory with raw SHARE-seq files')
    parser.add_argument('--output_dir', required=True, help='Directory to save processed data')
    parser.add_argument('--min_cells', type=int, default=2, help='Minimum cells per gene')
    parser.add_argument('--min_genes', type=int, default=100, help='Minimum genes per cell')
    args = parser.parse_args()

    # Define file paths
    RNA_COUNTS_FILE = args.input_dir+"/GSM4156610_brain.rna.counts.txt.gz"
    ATAC_COUNTS_FILE = args.input_dir+"/GSM4156599_brain.counts.txt.gz"
    ATAC_PEAKS_FILE = args.input_dir+"/GSM4156599_brain.peaks.bed.gz"
    CELLTYPE_FILE = args.input_dir+"/GSM4156599_brain_celltype.txt.gz"
    ATAC_BARCODES_FILE = args.input_dir+"/GSM4156599_brain.barcodes.txt.gz"
    ATAC_FRAGMENTS_FILE = args.input_dir+"/GSM4156599_brain.atac.fragments.sorted.bed.gz"
    ATAC_UNSORTED_FRAGMENTS_FILE = args.input_dir+"/GSM4156599_brain.atac.fragments.bed.gz"

    # Load data
    df_celltype = pd.read_csv(filepath_or_buffer=CELLTYPE_FILE, sep='\t')
    df_rna_counts = pd.read_csv(RNA_COUNTS_FILE, sep='\t', index_col=0)
    atac = sc.read_mtx(ATAC_COUNTS_FILE).T
    barcodes = pd.read_csv(ATAC_BARCODES_FILE, header=None)
    peaks = pd.read_csv(ATAC_PEAKS_FILE, sep="\t", header=None, names=["chrom", "start", "end"])

    # Set log file
    # Search for a directory "logs", if the directory does not exist, create it. The directory should be in the working dir
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / "prepare_shareseq_data.log"

    # Initialize adata_rna as empty ad.AnnData object
    adata_rna = ad.AnnData()

    # Process data
    adata_rna = get_rna_anndata(df_rna_counts, df_celltype)
    adata_rna = filter_rna_anndata(adata_rna, min_cells=args.min_cells, min_genes=args.min_genes, log_file=log_file)

    atac = process_atac_data(df_celltype, atac, barcodes, peaks)
    mu.atac.tl.locate_fragments(atac, str(ATAC_FRAGMENTS_FILE))
    pp_data = add_metrics_atac(atac, log_file=log_file)
    

    # Write data
    ad_rna_path = args.output_dir+"/c_sc_ad.h5ad"
    adata_rna.write_h5ad(ad_rna_path)
    atac_path = args.output_dir+"/c_atac_ad.h5ad"
    pp_data.write_h5ad(atac_path)
    