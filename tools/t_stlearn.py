from . import dataloader as dl
import stlearn as st
import argparse
from pathlib import Path

def preprocessing(adata):

    # Spatially-aware normalization and PCA
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]

    adata.obsm["X_morphology"] = adata.obsm["features"]

    st.em.run_pca(adata, n_comps=50)

    st.spatial.SME.SME_normalize(adata, use_data="raw")
    adata.X = adata.obsm['raw_SME_normalized']
    st.pp.scale(adata)
    st.em.run_pca(adata, n_comps=50)
    
    return adata


def cluster(adata, cluster_method, k_clusters=15, resolution=1.0):
    
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
    parser.add_argument("--cluster_method", type=str, default="louvain", help="Clustering method: kmeans, leiden, louvain.")
    parser.add_argument("--k_clusters", type=int, default=15, help="Number of clusters for k-means.")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for Louvain and Leiden clustering.")
    parser.add_argument("--do_normlog", type=bool, default=True, help="Whether to perform preprocessing.")
    
    args = parser.parse_args()

    adata = dl.open_anndata(args.input)
    
    adata = dl.shared_preprocessing(adata, args.do_normlog)

    adata = preprocessing(adata)
    adata = cluster(adata, cluster_method=args.cluster_method, k_clusters=args.k_clusters, resolution=args.resolution)

    dl.save_cluster(args.output, adata, "stlearn_cluster")
