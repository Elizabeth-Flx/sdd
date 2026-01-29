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

micromamba run -n "$multiomics_env" \
    python3 multiomics/writeX_tocsv.py --input_files "$rna_file" "$atac_file" --output_dir "$tmp_dir"

# 2. Run Seurat multi-omics integration in R
#micromamba run -n "$seurat_env" \
#    Rscript multiomics/seurat.R "$rna_csv" "$atac_csv" "$integrated_matrix_file"

# 3. Convert integrated matrix CSV back to AnnData
micromamba run -n "$multiomics_env" \
    python3 multiomics/build_integrated_anndata.py --seurat_matrix "multiomics/integrated_multiomics_matrix2.csv" --rna_h5ad "$rna_file" --atac_h5ad "$atac_file" --output_h5ad "$out_file"