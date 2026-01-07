import stlearn as st
import argparse
import scanpy as sc



def preprocessing(adata):
    # Preprocessing
    st.pp.filter_genes(adata, min_cells=3)
    st.pp.normalize_total(adata, target_sum=1e4)
    st.pp.log1p(adata)
    st.pp.scale(adata, max_value=10)

    # Spatially-aware normalization and PCA
    print(adata.obsm.keys())
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
    
    # adata.obsm["X_morphology"] = adata.obsm["features"]

    st.em.run_pca(adata, n_comps=50)

    # st.spatial.SME.SME_normalize(adata, use_data="raw")
    # adata.X = adata.obsm['raw_SME_normalized']
    # st.pp.scale(adata)
    # st.em.run_pca(adata, n_comps=50)
    
    return adata


def preprocessing_atac(adata):
    # Preprocessing for ATAC data
    st.pp.filter_genes(adata, min_cells=3)
    # st.pp.normalize_total(adata, target_sum=1e4)
    st.pp.scale(adata, max_value=10)

    # PCA
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
    st.em.run_pca(adata, n_comps=50)

    # st.spatial.SME.SME_normalize(adata, use_data="raw")
    # adata.X = adata.obsm['raw_SME_normalized']
    # st.pp.scale(adata)
    # st.em.run_pca(adata, n_comps=50)
    
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
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to input AnnData h5ad file")
    parser.add_argument("--output_h5ad", type=str, required=False, help="Path to save output AnnData h5ad file")
    parser.add_argument("--n_clusters", type=int, default=19, help="Number of clusters for k-means.")
    parser.add_argument("--louvain_resolution", type=float, default=1.19, help="Resolution parameter for Louvain clustering.")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_h5ad)

    if 'atac' in args.input_h5ad.lower():
        adata = preprocessing_atac(adata)
    else:
        adata = preprocessing(adata)
    adata = cluster(adata, n_clusters=args.n_clusters, louvain_resolution=args.louvain_resolution)

    output = args.output_h5ad # if args.output_h5ad else "../files/stlearn_" + args.input_h5ad

    adata.write_h5ad(output)