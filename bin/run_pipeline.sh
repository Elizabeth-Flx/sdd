#!/bin/bash

# Usage: ./run_pipeline.sh <input_dir> <output_dir> <graphst_env> <compare_env>
# Example: ./run_pipeline.sh data t_graphst sdd_tangram

DATA_DIR="$1"
GRAPHST_ENV="$2"
SCANIT_ENV="$3"
STLEARN_ENV="$4"
VIS_ENV="$5"
COMPARE_ENV="$6"

IN_DIR="$DATA_DIR/shareseq"
SHARESEQ_PROC_DIR="$DATA_DIR/shareseq_processed"
TANGRAM_OUT_DIR="$DATA_DIR/tangram_output"
TOOLS_OUT_DIR="$DATA_DIR/tools"
TOOL_NAMES=("graphst" "scanit" "stlearn")

# Tool name - cluster key mapping, one tool can have multiple cluster keys
declare -A TOOL_CLUSTER_KEYS
TOOL_CLUSTER_KEYS["graphst"]="mclust"
TOOL_CLUSTER_KEYS["scanit"]="scanit_louvain"
TOOL_CLUSTER_KEYS["stlearn"]="stlearn_kmeans"
TOOL_CLUSTER_KEYS["stlearn"]="stlearn_louvain"


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


# 2. Run Tools
# GraphST on spatial and Tangram-mapped data
echo "Running GraphST on spatial data..."
TOOL_NAME="graphst"
echo "Activating GraphST conda environment: $GRAPHST_ENV"
source ~/miniforge3/bin/activate "$GRAPHST_ENV"

python3 bin/t_graph_st.py --input_h5ad "$DATA_DIR/DLPFC_151671_processed.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_dlpfc_graphst.h5ad"
python3 bin/t_graph_st.py --input_h5ad "$TANGRAM_OUT_DIR/c_sp_ad.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_graphst.h5ad"
echo "Running GraphST on Tangram-mapped data..."
python3 bin/t_graph_st.py --input_h5ad "$TANGRAM_OUT_DIR/sp_atac.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_atac_graphst.h5ad"

# Scanit on on spatial and Tangram-mapped data
TOOL_NAME="scanit"
echo "Running Scanit on spatial data..."
echo "Activating Scanit conda environment: $SCANIT_ENV"
source ~/miniforge3/bin/activate "$SCANIT_ENV"

python3 bin/t_scanit.py --input_h5ad "$DATA_DIR/DLPFC_151671_processed.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_dlpfc_scanit.h5ad"
python3 bin/t_scanit.py --input_h5ad "$TANGRAM_OUT_DIR/c_sp_ad.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_scanit.h5ad"
echo "Running Scanit on Tangram-mapped data..."
python3 bin/t_scanit.py --input_h5ad "$TANGRAM_OUT_DIR/sp_atac.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_atac_scanit.h5ad"

# Stlearn on on spatial and Tangram-mapped data
TOOL_NAME="stlearn"
echo "Running Stlearn on spatial data..."
echo "Activating Stlearn conda environment: $STLEARN_ENV"
source ~/miniforge3/bin/activate "$STLEARN_ENV"
python3 bin/t_stlearn.py --input_h5ad "$DATA_DIR/DLPFC_151671_processed.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_dlpfc_stlearn.h5ad"
python3 bin/t_stlearn.py --input_h5ad "$TANGRAM_OUT_DIR/c_sp_ad.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_stlearn.h5ad"
echo "Running Stlearn on Tangram-mapped data..."
python3 bin/t_stlearn.py --input_h5ad "$TANGRAM_OUT_DIR/sp_atac.h5ad" --output_h5ad "$TOOLS_OUT_DIR/$TOOL_NAME/sp_atac_stlearn.h5ad"


# Switch to visualization environment
source ~/miniforge3/bin/activate "$VIS_ENV"

# 3. Visualizations
echo "Generating visualizations..."
# Iterate over mapping
for TOOL_NAME in "${TOOL_NAMES[@]}"; do
    # If the tool is stlearn, we have to switch back to the stlearn environment
    if [ "$TOOL_NAME" == "stlearn" ]; then
        echo "Activating Stlearn conda environment for visualization: $STLEARN_ENV"
        source ~/miniforge3/bin/activate "$STLEARN_ENV"
    fi
    # For each cluster key associated with the tool
    for key in ${TOOL_CLUSTER_KEYS[$TOOL_NAME]}; do
        echo "Generating visualizations for tool: $TOOL_NAME with cluster key: $key"
        python3 bin/vis_tools.py --input_dir "$TOOLS_OUT_DIR/$TOOL_NAME" --output_dir "$TOOLS_OUT_DIR/$TOOL_NAME/figures" --clusterkey "$key"
    done
    # Switch back to visualization environment if we switched to stlearn
    if [ "$TOOL_NAME" == "stlearn" ]; then
        echo "Switching back to visualization conda environment: $VIS_ENV"
        source ~/miniforge3/bin/activate "$VIS_ENV"
    fi
done

# 4. (Optional) Switch to comparison environment
# echo "Activating comparison conda environment: $COMPARE_ENV"
# source ~/miniforge3/bin/activate "$COMPARE_ENV"

# 5. Run comparison
# echo "Running clustering comparison..."
# for TOOL_NAME in "${TOOL_NAMES[@]}"; do
#     echo "Comparing results for tool: $TOOL_NAME"
#     python3 bin/compare_clustering.py --input_dir "$TOOLS_OUT_DIR/$TOOL_NAME" --output_dir "$TOOLS_OUT_DIR/$TOOL_NAME/compare_results" --tool "$TOOL_NAME"
# done

echo "Pipeline complete."