from . import dataloader as dl
import stlearn as st
import argparse
from pathlib import Path
import numpy as np

def preprocessing(adata):
    
    #library_id = list(adata.uns["spatial"].keys())[0]  # take the first library
    #adata.uns['spatial'][library_id]['use_quality'] = 'hires'  # or False
    
    #print(adata.uns["spatial"][library_id]["images"].keys())
    #print(adata.uns["spatial"][library_id].keys())
    
    #print(adata.obs.keys())
    
    try:
        adata.obs['imagerow'] = adata.obs['array_row']
        adata.obs['imagecol'] = adata.obs['array_col']
    except:
        pass

    try:
        adata.obs['imagerow'] = adata.obs['centroid_y']
        adata.obs['imagecol'] = adata.obs['centroid_x']
    except:
        pass

    try:
        adata.obs['array_row'] = adata.obs['centroid_y']
        adata.obs['array_col'] = adata.obs['centroid_x']
    except:
        pass
    # Spatially-aware normalization and PCA
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]

    # only add if not already present
    try:
        if "X_morphology" not in adata.obsm:
            adata.obsm["X_morphology"] = adata.obsm["features"]
    except:
        adata.obsm["X_morphology"] = np.zeros((adata.n_obs, 1))

    st.em.run_pca(adata, n_comps=50)

    st.spatial.SME.SME_normalize(adata, use_data="raw")
    adata.X = adata.obsm['raw_SME_normalized']
    st.pp.scale(adata)
    st.em.run_pca(adata, n_comps=50)
    
    return adata


def cluster(adata, cluster_method, k_clusters=12, resolution=1.0):
    
    if cluster_method == "kmeans":
        st.tl.clustering.kmeans(adata, n_clusters=k_clusters, use_data="X_pca", key_added="stlearn_cluster")
    elif cluster_method == "louvain":
        st.pp.neighbors(adata, n_neighbors=17, use_rep='X_pca')
        st.tl.clustering.louvain(adata, resolution=resolution, key_added="stlearn_cluster")
    elif cluster_method == "leiden":
        st.pp.neighbors(adata, n_neighbors=17, use_rep='X_pca')
        st.tl.clustering.leiden(adata, resolution=resolution, key_added="stlearn_cluster")

    return adata



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--output", type=str, required=True, help="Output csv file")    
    
    # parser.add_argument("--input", type=str, required=True, help="Name of the input AnnData file.")
    parser.add_argument("--cluster_method", type=str, default="kmeans", help="Clustering method: kmeans, leiden, louvain.")
    parser.add_argument("--k_clusters", type=int, default=12, help="Number of clusters for k-means.")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for Louvain and Leiden clustering.")
    parser.add_argument("--do_normlog", type=bool, default=True, help="Whether to perform preprocessing.")
    
    args = parser.parse_args()

    adata = dl.open_anndata(args.input)
    
    adata = dl.shared_preprocessing(adata, args.do_normlog)

    adata = preprocessing(adata)
    adata = cluster(adata, cluster_method=args.cluster_method, k_clusters=args.k_clusters, resolution=args.resolution)

    dl.save_cluster(args.output, adata, "stlearn_cluster")
