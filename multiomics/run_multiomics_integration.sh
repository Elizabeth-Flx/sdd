#! /bin/bash


# Function call: bash bin/run_multiomics_integration.sh rna_file atac_file multiomics_env seurat_env tmp_dir out_file

rna_file="$1"
atac_file="$2"
multiomics_env="$3"
seurat_env="$4"
tmp_dir="$5"
out_file="$6"

# The output of the first script is the filenames with _counts.csv suffix in tmp_dir
rna_csv="$tmp_dir/$(basename "${rna_file%.*}")_counts.csv"
atac_csv="$tmp_dir/$(basename "${atac_file%.*}")_counts.csv"
# Integrated matrix output file, temporary
integrated_matrix_file="$tmp_dir/integrated_multiomics_matrix.csv"

# Activate multi-omics integration environment
source ~/miniforge3/bin/activate "$multiomics_env"
# 1. Convert AnnData.X to CSV for Seurat
python3 bin/writeX_tocsv.py --input_files "$rna_file" "$atac_file" --output_dir "$tmp_dir"
# Deactivate multi-omics integration environment
conda deactivate
# 2. Run Seurat multi-omics integration in R
source ~/miniforge3/bin/activate "$seurat_env"
Rscript bin/seurat.R "$rna_csv" "$atac_csv" "$integrated_matrix_file"
# Deactivate Seurat environment
conda deactivate

# 3. Convert integrated matrix CSV back to AnnData
source ~/miniforge3/bin/activate "$multiomics_env"
python3 bin/build_integrated_anndata.py --seurat_matrix "$integrated_matrix_file" --rna_h5ad "$rna_file" --atac_h5ad "$atac_file" --output_h5ad "$out_file"