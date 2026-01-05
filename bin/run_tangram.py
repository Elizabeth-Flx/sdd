#! /usr/bin/env python3

import scanpy as sc
import pandas as pd
import tangram as tg
import argparse as ap
import matplotlib.pyplot as plt
import os
import shutil

NUM_EPOCHS = 500
DEVICE = "cpu"


def run_tangram(ad_rna, ad_sp, ad_atac_genes):
    # Ensure gene sets are identical
    genes_rna = set(ad_rna.var_names)
    genes_sp = set(ad_sp.var_names)
    # Adapt to common genes only
    common_genes = genes_rna.intersection(genes_sp)
    ad_rna = ad_rna[:, list(common_genes)].copy()
    ad_sp = ad_sp[:, list(common_genes)].copy()
    
    genes_rna = set(ad_rna.var_names)
    genes_sp = set(ad_sp.var_names)

    # cell_label = 'celltype'
    
    if genes_rna != genes_sp:
        raise ValueError("Gene sets in ad_rna and ad_sp are not identical. Please check preprocessing.")

    tg.pp_adatas(ad_rna, ad_sp, genes=ad_rna.var_names)

    ad_map = tg.map_cells_to_space(
        ad_rna, ad_sp,
        mode="cells",
        density_prior='rna_count_based',
        num_epochs=NUM_EPOCHS,
        device=DEVICE
    )

    ad_atac_genes.obs_names = ad_atac_genes.obs["barcode"].values
    ad_atac_genes.obs.index = ad_atac_genes.obs_names
    common_cells = ad_map.obs_names.intersection(ad_atac_genes.obs_names)
    common_cells = sorted(list(common_cells))  # Ensure same order
    shared_index = pd.Index(common_cells, dtype="str", name=None)
    print(f"Number of common cells for gene projection: {len(shared_index)}")
    ad_map_sub = ad_map[shared_index].copy()
    ad_atac_genes_sub = ad_atac_genes[shared_index].copy()

    ad_map_sub.obs_names = shared_index
    ad_atac_genes_sub.obs_names = shared_index
    ad_map_sub.obs.index = shared_index
    ad_atac_genes_sub.obs.index = shared_index
    
    ad_ge = tg.project_genes(ad_map_sub, ad_atac_genes_sub)

    tg.project_cell_annotations(ad_map_sub, ad_ge, "celltype")
    ad_ge.obsm['spatial'] = ad_sp.obsm['spatial']
    annotation_list = list(pd.unique(ad_rna.obs['celltype']))

    return ad_ge, annotation_list


# def gene_projection(ad_map, ad_atac_genes):
#     # Change obs of ad_atac_genes to the barcode column in obs
#     ad_atac_genes.obs_names = ad_atac_genes.obs["barcode"].values
#     ad_atac_genes.obs.index = ad_atac_genes.obs_names
#     common_cells = ad_map.obs_names.intersection(ad_atac_genes.obs_names)
#     common_cells = sorted(list(common_cells))  # Ensure same order
#     shared_index = pd.Index(common_cells, dtype="str", name=None)
#     print(f"Number of common cells for gene projection: {len(shared_index)}")
#     ad_map_sub = ad_map[shared_index].copy()
#     ad_atac_genes_sub = ad_atac_genes[shared_index].copy()

#     ad_map_sub.obs_names = shared_index
#     ad_atac_genes_sub.obs_names = shared_index
#     ad_map_sub.obs.index = shared_index
#     ad_atac_genes_sub.obs.index = shared_index
#     # Double-check
#     print("All index values match:", all(ad_map_sub.obs.index == ad_atac_genes_sub.obs.index))
#     print("Index objects are the same:", ad_map_sub.obs.index is ad_atac_genes_sub.obs.index)
#     print("obs_names objects are the same:", ad_map_sub.obs_names is ad_atac_genes_sub.obs_names)

#     tg.utils.project_genes(ad_map_sub, ad_atac_genes_sub)
#     return ad_map_sub


def plot_tangram_results(ad_map, ad_sp, ad_rna, annotation_list, plot_path):

    tg.plot_cell_annotation_sc(ad_map, annotation_list, x='x', y='y',spot_size= 60, scale_factor=0.1, perc=0.02)   # type:ignore
    plt.tight_layout()
    plt.savefig(plot_path+"celltype_spatial_annotation.png", dpi=300)
    
    tg.plot_training_scores(ad_map, bins=20, alpha=.5)
    plt.tight_layout()
    plt.savefig(plot_path+"tangram_training_scores.png", dpi=300)
    
    anno = list(pd.unique(ad_rna.obs['celltype']))
    tg.plot_cell_annotation_sc(ad_sp, anno, spot_size=40, perc=0.02)  # type:ignore
    plt.tight_layout()
    plt.savefig(plot_path+"tangram_celltype_spatial_annotation.png", dpi=300)
    

if __name__ == "__main__":

    # Function call: python3 bin/run_tangram.py --input_dir data/shareseq_processed --output_dir data/tangram_output

    parser = ap.ArgumentParser(description="Run Tangram mapping for SHARE-seq data.")
    parser.add_argument('--input_dir', required=True, help='Directory with processed AnnData files')
    parser.add_argument('--output_dir', required=True, help='Directory to save Tangram results')
    args = parser.parse_args()

    ad_sp = sc.read_h5ad(args.input_dir+"/c_sp_ad.h5ad")
    ad_atac = sc.read_h5ad(args.input_dir+"/c_atac_ad.h5ad")
    ad_rna = sc.read_h5ad(args.input_dir+"/c_sc_ad.h5ad")
    ad_atac_genes = sc.read_h5ad(args.input_dir+"/c_atac_genes_ad.h5ad")

    # Run Tangram mapping
    ad_map_ge, annotation_list = run_tangram(ad_rna, ad_sp, ad_atac_genes)

    # Save ad_map
    output_path = f"{args.output_dir}/sp_atac.h5ad"
    ad_map_ge.write_h5ad(output_path)  # type:ignore

    # Copy sp file to output dir
    output_spatial_path = f"{args.output_dir}/c_sp_ad.h5ad"
    shutil.copy(f"{args.input_dir}/c_sp_ad.h5ad", output_spatial_path)

    plot_path = f"{args.output_dir}/figures/"
    # Create directory if it doesn't exist
    os.makedirs(plot_path, exist_ok=True)
    # plot_tangram_results(ad_map, ad_sp, ad_rna, annotation_list, plot_path)