#! /usr/bin/env python3

import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
import muon as mu
from pathlib import Path
import sys
import re
import argparse
import genetic_annotation as ga
import shutil


## Set constants
MIN_CELLS = 2
MIN_GENES = 100
N_TOP_GENES = 2000
FLAVOR = "seurat_v3"
FRAGMENT_MIN = 165
FRAGMENT_MAX = 12000
NUCLEOSOME_SIGNAL_THRESHOLD = 2
NUCLEOSOME_SIGNAL = 2000000

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


def filter_rna_anndata(adata_rna, log_file=None):

    mito_count = adata_rna.var_names.str.startswith("Mt-").sum()
    # Write mito count to log file
    if log_file:
        with open(log_file, "a") as f:
            f.write("FILTERING RNA ANNDATA\n")
            f.write(f"Mitochondrial genes count: {mito_count}\n")
            # Write parameters used for filtering
            f.write(f"min_cells: {MIN_CELLS}\n\
                        min_genes: {MIN_GENES}\n")

    sc.pp.calculate_qc_metrics(
        adata_rna, inplace=True, log1p=True
    )

    # Filter data
    sc.pp.filter_cells(adata_rna, min_genes=MIN_GENES) # pyright: ignore[reportCallIssue]
    sc.pp.filter_genes(adata_rna, min_cells=MIN_CELLS) # pyright: ignore[reportCallIssue]

    # Saving count data
    adata_rna.layers["counts"] = adata_rna.X.copy()
    sc.pp.normalize_total(adata_rna) # pyright: ignore[reportCallIssue]
    sc.pp.log1p(adata_rna)
    sc.pp.highly_variable_genes(adata_rna, n_top_genes=N_TOP_GENES, flavor=FLAVOR) # pyright: ignore[reportCallIssue]
    sc.tl.pca(adata_rna)

    return adata_rna


def process_atac_data(celltype, atac, barcodes, peaks):
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


def add_metrics_atac(atac, log_file=None):
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
    mu.atac.tl.nucleosome_signal(atac, n=NUCLEOSOME_SIGNAL)

    # Add group labels for above and below the nucleosome signal threshold
    atac.obs["nuc_signal_filter"] = [
        "NS_FAIL" if ns > NUCLEOSOME_SIGNAL_THRESHOLD else "NS_PASS"
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
        lambda x:  (x >= FRAGMENT_MIN) &(x <= FRAGMENT_MAX),  # parameter for filtering
    )
    mu.pp.filter_var(pp_data, "n_cells_by_counts", lambda x: x >= MIN_CELLS)
    if log_file is not None:
        with open(log_file, "a") as f:
            f.write(f"Number of cells after filtering on total_fragment_counts & n_cells_by_counts: {pp_data.n_obs}\n")
            
    return pp_data


def get_atac_genes(adata_rna, atac, input_dir):
    # Define the GTF URL
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz"

    # Load and process the genes, excluding chrY
    try:
        gene_coords_gencode = ga.load_gencode_gtf_genes(gtf_url, feature_type='gene', exclude_chroms=['chrY'])
    except Exception as e:
        print(f"Failed to load or process GTF: {e}")

    gene_coords_gencode = gene_coords_gencode[~gene_coords_gencode["chrom"].str.contains("chrY")]
    common_genes_gencode = adata_rna.var_names.intersection(gene_coords_gencode.index)
    adata_rna_annotated = adata_rna[:, list(common_genes_gencode)].copy()
        
    gene_coords_gencode_reset = gene_coords_gencode.reset_index()

    gene_coords_gencode_filtered = (
        gene_coords_gencode_reset
        .loc[gene_coords_gencode_reset['gene_name'].isin(common_genes_gencode)]
        .groupby(['gene_name', 'chrom'])
        .agg({
            'start': 'min',
            'end': 'max',
            'Strand': 'first',
            'gene_ids': 'first'
        })
        .reset_index()
    )

    gene_coords_gencode_filtered = gene_coords_gencode_filtered[
        ~((gene_coords_gencode_filtered['gene_name'] == 'Fth1') &
        (gene_coords_gencode_filtered['chrom'] == 'chr8'))
    ]

    gene_coords_gencode_filtered = gene_coords_gencode_filtered.set_index('gene_name')

    mdata = mu.MuData({'rna': adata_rna_annotated, 'atac': atac})
    mu.atac.tl.initialise_default_files(mdata, path=input_dir)
    adata_atac_genes = mu.atac.tl.count_fragments_features(
        mdata,
        features=gene_coords_gencode_filtered,
        stranded=True,
        count_reads=False
    )

    sc.pp.filter_genes(adata_atac_genes, min_cells=MIN_CELLS) # type: ignore
    sc.pp.normalize_total(adata_atac_genes) # type: ignore
    # Logarithmize the data
    sc.pp.log1p(adata_atac_genes)

    return adata_atac_genes


    
if __name__ == "__main__":

    # Function call: python3 bin/prepare_shareseq_data.py --input_dir data/shareseq --output_dir data/shareseq_processed

    parser = argparse.ArgumentParser(description="Prepare SHARE-seq data for downstream analysis.")
    parser.add_argument('--input_dir', required=True, help='Directory with raw SHARE-seq files')
    parser.add_argument('--output_dir', required=True, help='Directory to save processed data')
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
    atac_mtx = sc.read_mtx(ATAC_COUNTS_FILE).T
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
    adata_rna = filter_rna_anndata(adata_rna, log_file=log_file)

    atac_processed = process_atac_data(df_celltype, atac_mtx, barcodes, peaks)
    mu.atac.tl.locate_fragments(atac_processed, str(ATAC_FRAGMENTS_FILE))
    pp_data = add_metrics_atac(atac_processed, log_file=log_file)
    atac_genes = get_atac_genes(adata_rna, pp_data, args.input_dir)


    # Write processed RNA and ATAC data
    ad_rna_path = args.output_dir+"/c_sc_ad.h5ad"
    adata_rna.write_h5ad(ad_rna_path)
    atac_path = args.output_dir+"/c_atac_ad.h5ad"
    pp_data.write_h5ad(atac_path)
    atac_gene_path = args.output_dir+"/c_atac_genes_ad.h5ad"
    atac_genes.write_h5ad(atac_gene_path)

    # Process and save spatial data subset to HVGs (if spatial file exists)
    spatial_path = Path(args.input_dir) / "c_sp_ad.h5ad"
    # Copy file to output dir
    output_spatial_path = Path(args.output_dir) / "c_sp_ad.h5ad"
    shutil.copy(spatial_path, output_spatial_path)
    