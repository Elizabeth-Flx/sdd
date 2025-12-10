#! /usr/bin/env python3

# NOTE: This is still a dummy script, changes will likely be necessary (exact formats etc. are still unknown)

import numpy as np
import argparse as ap
import scanpy as sc


"""
Each tool performs spatial domain detection on both spatial and spatial+ATAC data,
producing two AnnData objects. This script takes both AnnData objects and computes
a consensus annotation for each spatial spot based on the annotations from both objects.
"""

def get_consensus_annotation(sp, sp_atac):

    """
    Get consensus annotations between spatial and spatial+ATAC data.
    @param sp: AnnData object for spatial data
    @param sp_atac: AnnData object for spatial+ATAC data
    @return: AnnData object with consensus annotations
    """
    # Ensure both AnnData objects have the same spots
    if not np.array_equal(sp.obs_names, sp_atac.obs_names):
        raise ValueError("The spatial spots in both AnnData objects do not match.")
    consensus = sc.AnnData(sp.X, obs=sp.obs.copy(), var=sp.var.copy(), uns=sp.uns.copy())
    consensus.obs['consensus_annotation'] = np.where(
        sp.obs['celltype'] == sp_atac.obs['celltype'],
        sp.obs['celltype'],
        'ambiguous'
    )
    
    
    return consensus



if __name__ == "__main__":
    
    parser = ap.ArgumentParser(description="Find consensus annotations between spatial and spatial+ATAC data.")
    parser.add_argument('--sp_file', help='Path to the spatial AnnData file')
    parser.add_argument('--sp_atac_file', help='Path to the spatial+ATAC AnnData file')
    parser.add_argument('--output_file', help='Output file to save the consensus annotations')
    args = parser.parse_args()

    # Load AnnData objects
    ad_sp = sc.read_h5ad(args.sp_file)
    ad_sp_atac = sc.read_h5ad(args.sp_atac_file)

    # Get consensus annotations
    consensus_annotations = get_consensus_annotation(ad_sp, ad_sp_atac)

    # Save consensus annotations to file as h5ad
    consensus_annotations.write_h5ad(args.output_file)