#! /usr/bin/env python3

import re
import scanpy as sc
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import tangram as tg
import muon as mu
import anndata as ad

import os
import argparse as ap


def run_tangram(ad_rna, ad_sp, num_epochs=500, device='cpu'):
    overlapping_genes: set[str] = set(ad_rna.var_names.str.lower()) & set(ad_sp.var_names.str.lower())

    tg.pp_adatas(ad_rna, ad_sp, genes=overlapping_genes)

    ad_map = tg.map_cells_to_space(ad_rna, ad_sp,
        mode="cells",
        density_prior='rna_count_based',
        num_epochs=num_epochs,
        # device="cuda:0",
        device=device,
    )
    tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype")
    annotation_list = list(pd.unique(ad_rna.obs['celltype']))

    return ad_map, annotation_list
    


if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Run Tangram mapping for SHARE-seq data.")
    parser.add_argument('--input_dir', required=True, help='Directory with processed AnnData files')
    parser.add_argument('--output_dir', required=True, help='Directory to save Tangram results')
    parser.add_argument('--tool_name', required=True, help='Spatial domain detection tool name')
    parser.add_argument('--num_epochs', type=int, default=500, help='Number of epochs for Tangram')
    parser.add_argument('--device', default='cpu', help='Device for Tangram (cpu/cuda)')
    args = parser.parse_args()

    ad_sp = sc.read_h5ad(args.input_dir+"/c_sp_ad.h5ad")
    ad_rna = sc.read_h5ad(args.input_dir+"/c_sc_ad.h5ad")
    # ad_atac_genes = sc.read_h5ad(args.input_dir+"/c_atac_genes_ad.h5ad")

    # Run Tangram mapping
    ad_map, annotation_list = run_tangram(ad_rna, ad_sp, num_epochs=args.num_epochs, device=args.device)

    # TODO: Continue here