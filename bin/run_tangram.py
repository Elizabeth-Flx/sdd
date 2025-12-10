#! /usr/bin/env python3

import scanpy as sc
import pandas as pd
import tangram as tg
import argparse as ap

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

    tg.utils.project_genes(ad_map, ad_atac_genes, annotation_list)

    # Save ad_map
    output_path = f"{args.output_dir}/sp_atac.h5ad"
    ad_map.write_h5ad(output_path)