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


def cluster(adata, n_clusters=19, louvain_resolution=1.19):
    # Clustering with k-means
    st.tl.clustering.kmeans(adata, n_clusters=n_clusters, use_data="X_pca", key_added="stlearn_kmeans")

    # Clustering with Louvain
    st.pp.neighbors(adata, n_neighbors=17, use_rep='X_pca')
    st.tl.clustering.louvain(adata, resolution=louvain_resolution, key_added="stlearn_louvain")

    return adata



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell spatial transcriptomics data with STLearn.")
    
    parser.add_argument("--input", type=str, required=True, help="Input anndata file")
    parser.add_argument("--output", type=str, required=True, help="Output csv file")    
    
    # parser.add_argument("--input", type=str, required=True, help="Name of the input AnnData file.")
    parser.add_argument("--n_clusters", type=int, default=19, help="Number of clusters for k-means.")
    parser.add_argument("--louvain_resolution", type=float, default=1.19, help="Resolution parameter for Louvain clustering.")
    parser.add_argument("--do_normlog", type=bool, default=True, help="Whether to perform preprocessing.")
    
    args = parser.parse_args()

    adata = dl.open_anndata(args.input)
    
    adata = dl.shared_preprocessing(adata, args.do_normlog)

    adata = preprocessing(adata)
    adata = cluster(adata, n_clusters=args.n_clusters, louvain_resolution=args.louvain_resolution)

    dl.save_cluster(args.output, adata, "stlearn_louvain")
