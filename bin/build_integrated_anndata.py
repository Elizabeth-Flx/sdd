import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import os

def main(args):
    # Load original data
    sp_rna = sc.read_h5ad(args.rna_h5ad)
    sp_atac = sc.read_h5ad(args.atac_h5ad)
    seurat_matrix = pd.read_csv(args.seurat_matrix, index_col=0)

    # Print shapes
    # print('sp_rna_counts:', sp_rna_counts.shape)
    # print('sp_atac_counts:', sp_atac_counts.shape)
    # print('seurat_matrix:', seurat_matrix.shape)

    # Prepare obs
    obs_rna = sp_rna.obs.copy()
    obs_atac = sp_atac.obs.copy()
    obs_rna.index = [f"{idx}_1" for idx in obs_rna.index]
    obs_atac.index = [f"{idx}_2" for idx in obs_atac.index]
    obs_combined = pd.DataFrame()
    obs_combined = pd.concat([obs_rna, obs_atac], axis=0)
    obs_combined = obs_combined.loc[seurat_matrix.index]

    # Prepare obsm (spatial)
    obsm_combined = {}
    if 'spatial' in sp_rna.obsm and 'spatial' in sp_atac.obsm:
        spatial_rna = pd.DataFrame(sp_rna.obsm['spatial'], index=[f"{name}_1" for name in sp_rna.obs_names])
        spatial_atac = pd.DataFrame(sp_atac.obsm['spatial'], index=[f"{name}_2" for name in sp_atac.obs_names])
        spatial_combined = pd.concat([spatial_rna, spatial_atac], axis=0)
        spatial_combined = spatial_combined.loc[seurat_matrix.index]
        obsm_combined['spatial'] = spatial_combined.values

    # Prepare var
    var_rna = pd.DataFrame(index=sp_rna.var_names)
    var_atac = pd.DataFrame(index=sp_atac.var_names)
    common_vars = var_rna.index.intersection(var_atac.index)
    var_combined = pd.DataFrame(index=common_vars)
    all_cols = set(sp_rna.var.keys()).union(set(sp_atac.var.keys()))
    for col in all_cols:
        if col in sp_rna.var.keys() and col in sp_atac.var.keys():
            var_combined[f"rna_{col}"] = sp_rna.var.loc[common_vars, col]
            var_combined[f"atac_{col}"] = sp_atac.var.loc[common_vars, col]
        else:
            if col in sp_rna.var.keys():
                var_combined[f"rna_{col}"] = sp_rna.var.loc[common_vars, col]
            if col in sp_atac.var.keys():
                var_combined[f"atac_{col}"] = sp_atac.var.loc[common_vars, col]

    # Build AnnData object
    adata_seurat = sc.AnnData(
        X=seurat_matrix.values,
        obs=obs_combined,  # type: ignore
        var=var_combined,
        obsm=obsm_combined
    )
    print(adata_seurat)
    adata_seurat.write(args.output_h5ad)
    print(f"Saved integrated AnnData to {args.output_h5ad}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build integrated AnnData from Seurat and spatial data.")
    parser.add_argument('--rna_h5ad', type=str, required=True, help='Path to spatial RNA AnnData (.h5ad)')
    parser.add_argument('--atac_h5ad', type=str, required=True, help='Path to spatial ATAC AnnData (.h5ad)')
    parser.add_argument('--seurat_matrix', type=str, required=True, help='Path to integrated Seurat matrix CSV')
    parser.add_argument('--output_h5ad', type=str, required=True, help='Output path for integrated AnnData (.h5ad)')
    args = parser.parse_args()
    main(args)
