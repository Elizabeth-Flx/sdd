import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import *
from plotnine import *
import scanpy as sc
import anndata as ad
from matplotlib.path import Path
import os
from scib.metrics import *
from sklearn.neighbors import NearestNeighbors
from pathlib import Path


# loading csvs and proof of consistency across barcodes

def load_csvs_to_arrays(data_path, files, label_col="spatial_cluster", barcode_col="barcode"):
    dict_all_methods = {}
    groundtruth = None

    ref_barcodes = None  # barcode order we will enforce
    seen = 0

    for csv_file in sorted(files):
        if not csv_file.endswith(".csv"):
            continue

        name = os.path.splitext(csv_file)[0]
        fpath = os.path.join(data_path, csv_file)

        df = pd.read_csv(fpath)

        # basic validation
        if barcode_col not in df.columns or label_col not in df.columns:
            raise ValueError(
                f"{csv_file} must have columns '{barcode_col}' and '{label_col}'. "
                f"Found: {list(df.columns)}"
            )

        barcodes = df[barcode_col].astype(str).to_numpy()
        labels = df[label_col].astype(int).to_numpy()


        # check duplicates
        if len(np.unique(barcodes)) != len(barcodes):
            dup = pd.Series(barcodes).duplicated()
            examples = pd.Series(barcodes)[dup].unique()[:10]
            raise ValueError(f"{csv_file}: duplicated barcodes (e.g. {examples})")

        # enforce identical barcode order
        if ref_barcodes is None:
            ref_barcodes = barcodes
        else:
            if len(barcodes) != len(ref_barcodes):
                raise ValueError(f"{csv_file}: different number of rows than reference CSV")

            if not np.array_equal(barcodes, ref_barcodes):
                mismatch = np.where(barcodes != ref_barcodes)[0][0]
                raise ValueError(
                    f"{csv_file}: barcode order mismatch at row {mismatch}. "
                    f"Expected {ref_barcodes[mismatch]} but got {barcodes[mismatch]}."
                )

        if name.lower() == "groundtruth":
            groundtruth = labels
        else:
            dict_all_methods[name] = labels
        print(f"Loaded {name} from {csv_file} (n={len(labels)})")
        seen += 1
    if seen == 0:
        raise ValueError(f"No CSV files found in {data_path}")

    return dict_all_methods, groundtruth, ref_barcodes
        

# accuracy metrics helper functions

def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = pd.crosstab(y_true, y_pred)
    return np.sum(np.amax(contingency_matrix.values, axis=0)) / np.sum(contingency_matrix.values)


# accuracy: ARI, Purity, AMI

def calculate_accuracy_metrics(dict_all_methods, groundtruth = None):   
    if groundtruth is not None:
        # calculate ARI and AMI for each method in dict_all_methods
        ARI_scores = {}
        Purity_scores = {}
        AMI_scores = {}
        for method_name, values in dict_all_methods.items():
            ARI = adjusted_rand_score(groundtruth, values)
            Purity = purity_score(groundtruth, values)
            AMI = adjusted_mutual_info_score(groundtruth, values)
            ARI_scores[method_name] = ARI
            Purity_scores[method_name] = Purity
            AMI_scores[method_name] = AMI
        
        print(ARI_scores)
        print(Purity_scores)
        print(AMI_scores)
        return ARI_scores, Purity_scores, AMI_scores
    else:
        print("Groundtruth labels not found. Skipping accuracy metrics calculation.")
        

# spatial metrics helper functions: PAS and iLISI

def calculate_pas(spatial_coords, labels, k=10, threshold=0.6):
    n_spots = len(labels)
    
    # Find k nearest neighbors (including self)
    nn = NearestNeighbors(n_neighbors=k+1)
    nn.fit(spatial_coords)
    distances, indices = nn.kneighbors(spatial_coords)
    
    # Remove self from neighbors (first column)
    neighbor_indices = indices[:, 1:]
    ambiguous_count = 0
    
    for i in range(n_spots):
        # Get labels of k nearest neighbors
        neighbor_labels = labels[neighbor_indices[i]]
        current_label = labels[i]
        
        # Count neighbors with different label
        different_count = np.sum(neighbor_labels != current_label)
        
        # Check if ambiguous (≥60% different)
        if different_count >= k * threshold:
            ambiguous_count += 1
    
    return ambiguous_count / n_spots

# need to use spatial coordinates for spatial data
# assume it is saved in adata.obsm['spatial']

def calculate_spatial_metrics(adata, dict_all_methods):
    
    ilisi_scores = {}
    pas_scores = {}
    
    for method, labels in dict_all_methods.items():
        adata_temp = adata.copy()
        adata_temp.obs[method] = pd.Categorical(labels)

    # create knn graph based on spatial coordinates
        sc.pp.neighbors(adata_temp, use_rep='spatial', n_neighbors=90)
        
        ilisi_score = ilisi_graph(
            adata_temp,
            batch_key=method,  # Using method labels as "batches"
            type_='knn',     # Using nearest neighbor graph
            k0=90,            # Number of nearest neighbors
            scale=True        # Scale between 0 and 1
        )
        
        pas_score = calculate_pas(
            adata_temp.obsm['spatial'], 
            dict_all_methods[method], 
            k=10, 
            threshold=0.6)
        
        ilisi_scores[method] = ilisi_score
        pas_scores[method] = pas_score
        
        print(f"{method} iLISI score: {ilisi_score}")
        print(f"{method} PAS score: {pas_score}")
        
    return ilisi_scores, pas_scores


# silhouette score functions helper: intracluster and intercluster

# for transcriptome-based silhouette score
# Method 1: Using the raw expression data (adata.X)
# calculates in expression space, higher = better separation
def calculate_silhouette_expression(adata, method_name):
    labels = adata.obs[method_name].values
    
    # Use PCA embeddings for faster computation
    if 'X_pca' in adata.obsm:
        X = adata.obsm['X_pca'][:, :50]
    else:
        # Use raw expression
        X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    score = silhouette_score(
        X, 
        labels, 
        metric='euclidean',
        #sample_size=min(1000, len(labels))  # Subsample for speed
    )
    return score

# for physical space-based silhouette score
def calculate_silhouette_physical(adata, method_name):
    labels = adata.obs[method_name].values
    spatial_coords = adata.obsm['spatial']
    
    score = silhouette_score(
        spatial_coords, 
        labels, 
        metric='euclidean',
        #sample_size=min(1000, len(labels))  # Subsample for speed
    )
    return score

# both silhouette scores together

def calculate_silhouette_metrics(adata, dict_all_methods):
    silhouette_expression_scores = {}
    silhouette_physical_scores = {}
    
    for method in dict_all_methods.keys():
        expr_score = calculate_silhouette_expression(adata, method)
        phys_score = calculate_silhouette_physical(adata, method)
        
        silhouette_expression_scores[method] = expr_score
        silhouette_physical_scores[method] = phys_score
        
        print(f"{method} Silhouette Expression score: {expr_score}")
        print(f"{method} Silhouette Physical score: {phys_score}")
        
    return silhouette_expression_scores, silhouette_physical_scores

def generate_metrics_df_from_dicts(metric_dicts):
    datasets = sorted(set().union(*[set(d.keys()) for d in metric_dicts.values()]))
    data = []
    for dataset in datasets:
        # union of methods across all metrics for this dataset
        methods = set()
        for md in metric_dicts.values():
            methods |= set(md.get(dataset, {}).keys())
        methods = sorted(methods)

        for method in methods:
            row = {"Dataset": dataset, "Method": method}

            # fill every metric column (NaN if missing)
            for metric_name, md in metric_dicts.items():
                row[metric_name] = md.get(dataset, {}).get(method, np.nan)

            data.append(row)

    df_metrics = pd.DataFrame(data).sort_values(["Dataset", "Method"]).reset_index(drop=True)

    # generate mean dataset metrics
    metric_cols = list(metric_dicts.keys())  # ["ARI","Purity",...]
    mean_df_metrics = (
        df_metrics
        .groupby("Method", as_index=False)[metric_cols]
        .mean(numeric_only=True)
    )
    mean_df_metrics
    
    return df_metrics, mean_df_metrics