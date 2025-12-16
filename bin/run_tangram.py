#! /usr/bin/env python3

import scanpy as sc
import pandas as pd
import tangram as tg
import argparse as ap
import matplotlib.pyplot as plt
import os

NUM_EPOCHS = 500
DEVICE = "cpu"


def run_tangram(ad_rna, ad_sp):
    overlapping_genes: set[str] = set(ad_rna.var_names.str.lower()) & set(ad_sp.var_names.str.lower())

    tg.pp_adatas(ad_rna, ad_sp, genes=overlapping_genes)

    ad_map = tg.map_cells_to_space(ad_rna, ad_sp,
        mode="cells",
        density_prior='rna_count_based',
        num_epochs=NUM_EPOCHS,
        device=DEVICE,
    )
    tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype")
    annotation_list = list(pd.unique(ad_rna.obs['celltype']))

    return ad_map, annotation_list


def plot_tangram_results(ad_map, ad_sp, ad_rna, annotation_list, plot_path):

    tg.plot_cell_annotation_sc(ad_sp, annotation_list, x='x', y='y',spot_size= 60, scale_factor=0.1, perc=0.02)
    plt.tight_layout()
    plt.savefig(plot_path+"celltype_spatial_annotation.png", dpi=300)
    
    tg.plot_training_scores(ad_map, bins=20, alpha=.5)
    plt.tight_layout()
    plt.savefig(plot_path+"tangram_training_scores.png", dpi=300)
    
    ad_sp.obs["tangram_celltype"] = ad_sp.obsm["tangram_ct_pred"].idxmax(axis=1)
    anno = list(pd.unique(ad_rna.obs['celltype']))
    tg.plot_cell_annotation_sc(ad_sp, anno, spot_size=40, perc=0.02)
    plt.tight_layout()
    plt.savefig(plot_path+"tangram_celltype_spatial_annotation.png", dpi=300)
    

if __name__ == "__main__":

    # Function call: python3 bin/run_tangram.py --input_dir data/shareseq_processed --output_dir data/tangram_output

    parser = ap.ArgumentParser(description="Run Tangram mapping for SHARE-seq data.")
    parser.add_argument('--input_dir', required=True, help='Directory with processed AnnData files')
    parser.add_argument('--output_dir', required=True, help='Directory to save Tangram results')
    args = parser.parse_args()

    ad_sp = sc.read_h5ad(args.input_dir+"/c_sp_ad.h5ad")
    ad_rna = sc.read_h5ad(args.input_dir+"/c_sc_ad.h5ad")
    ad_atac_genes = sc.read_h5ad(args.input_dir+"/c_atac_genes_ad.h5ad")

    # Run Tangram mapping
    ad_map, annotation_list = run_tangram(ad_rna, ad_sp)

    # tg.utils.project_genes(ad_map, ad_atac_genes, annotation_list)

    # Save ad_map
    output_path = f"{args.output_dir}/sp_atac.h5ad"
    ad_map.write_h5ad(output_path)

    plot_path = f"{args.output_dir}/figures/"
    # Create directory if it doesn't exist
    os.makedirs(plot_path, exist_ok=True)
    plot_tangram_results(ad_map, ad_sp, ad_rna, annotation_list, plot_path)