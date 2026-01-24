#! /bin/bash

# Function call: bash bin/run_tools.sh data t_graphst t_scanit t_stlearn sdd_tangram

DATA_DIR="$1"
GRAPHST_ENV="$2"
SCANIT_ENV="$3"
STLEARN_ENV="$4"
VIS_ENV="$5"

TOOLS_OUT_DIR="$DATA_DIR/tools"
TOOL_NAMES=("graphst") #"graphst" "scanit" 

# matrix_all="$DATA_DIR/sp_integrated_all.h5ad"
matrix_rna="$DATA_DIR/sp_integrated_rna.h5ad"

# Run all tools for both datasets
for tool in "${TOOL_NAMES[@]}"; do
    echo "Running $tool on all modalities..."
    # Activate the corresponding environment
    case $tool in
        "graphst")
            source ~/miniforge3/bin/activate "$GRAPHST_ENV"
            ;;
        "scanit")
            source ~/miniforge3/bin/activate "$SCANIT_ENV"
            ;;
        "stlearn")
            source ~/miniforge3/bin/activate "$STLEARN_ENV"
            ;;
        *)
            echo "Unknown tool: $tool"
            exit 1
            ;;
    esac
    # Run the tool
    # python3 "bin/t_${tool}.py" --input_h5ad "$matrix_all" --output_h5ad "$TOOLS_OUT_DIR/${tool}/sp_integrated_all.h5ad"
    python3 "bin/t_${tool}.py" --input_h5ad "$matrix_rna" --output_h5ad "$TOOLS_OUT_DIR/${tool}/sp_integrated_rna.h5ad"
    # Deactivate the environment
    conda deactivate
    echo "$tool completed on all modalities."
done