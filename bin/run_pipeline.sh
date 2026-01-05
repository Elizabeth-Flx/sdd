#!/bin/bash

# Usage: ./run_pipeline.sh <input_dir> <output_dir> <graphst_env> <compare_env>
# Example: ./run_pipeline.sh data t_graphst sdd_tangram

DATA_DIR="$1"
GRAPHST_ENV="$2"
COMPARE_ENV="$3"

IN_DIR="$DATA_DIR/shareseq"
SHARESEQ_PROC_DIR="$DATA_DIR/shareseq_processed"
TANGRAM_OUT_DIR="$DATA_DIR/tangram_output"
TOOLS_OUT_DIR="$DATA_DIR/tools"

echo "Input Directory: $IN_DIR"
echo "Shareseq Processed Directory: $SHARESEQ_PROC_DIR"
echo "Tangram Output Directory: $TANGRAM_OUT_DIR"
echo "Tools Output Directory: $TOOLS_OUT_DIR"

# 1. Preprocessing (prepare_shareseq_data.py)
echo "Running preprocessing..."
python3 bin/prepare_shareseq_data.py --input_dir "$IN_DIR" --output_dir "$SHARESEQ_PROC_DIR"

echo "Run tangram mapping..."
# Run Tangram mapping
python3 bin/run_tangram.py --input_dir "$SHARESEQ_PROC_DIR" --output_dir "$TANGRAM_OUT_DIR"
# 2. (Optional) Activate GraphST environment
echo "Activating GraphST conda environment: $GRAPHST_ENV"
source ~/miniforge3/bin/activate "$GRAPHST_ENV"

# 3. Run GraphST on spatial and Tangram-mapped data
echo "Running GraphST on spatial data..."
TOOL_NAME="graphst"
python3 bin/t_graph_st.py --input_h5ad "$SHARESEQ_PROC_DIR/c_sp_ad.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/c_sp_ad_graphst.h5ad"
echo "Running GraphST on Tangram-mapped data..."
python3 bin/t_graph_st.py --input_h5ad "$TANGRAM_OUT_DIR/sp_atac.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_atac_graphst.h5ad"

# 4. (Optional) Switch to comparison environment
echo "Activating comparison conda environment: $COMPARE_ENV"
source ~/miniforge3/bin/activate "$COMPARE_ENV"

# 5. Run comparison
echo "Running clustering comparison..."
python3 bin/compare_clustering.py --input_dir "$TOOLS_OUT_DIR/$TOOL_NAME" --output_dir "$TOOLS_OUT_DIR/$TOOL_NAME/compare_results"

echo "Pipeline complete."