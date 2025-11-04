#!/bin/bash
#
# SLURM wrapper for Phase 6: Normalization and Covariate Adjustment
# This script executes the Python script to perform Inverse Normal Transformation (INT)
# on the edQTL features and prepares the final covariate matrix for FastQTL.

# --- SLURM Directives ---
# Set the job name
#SBATCH --job-name=edQTL_P6_Normalize
# Set the output file path (%j is the job ID)
#SBATCH --output=logs/phase6_processing_%j.out
# Set the error file path
#SBATCH --error=logs/phase6_processing_%j.err
# Set the partition (adjust this based on your cluster)
#SBATCH --partition=bigmem
# Reserve sufficient memory (assuming 6000 individuals * N_Features is large)
#SBATCH --mem=128G 
# Set time limit (adjust as needed)
#SBATCH --time=04:00:00

# --- Configuration ---
# Path to the Python script
SCRIPT="./normalize_and_covariate_phase6.py"

# Input file from Phase 5 (edQTL Feature Matrix - Raw ERs)
INPUT_FILE="./phase5_edQTL_features/edQTL_features_all_genes_all_celltypes_p5.tsv"

# Output directory for Phase 6
OUTPUT_DIR="./phase6_normalized_edQTL"

# Final output files for FastQTL
OUTPUT_NORMALIZED_FILE="${OUTPUT_DIR}/normalized_edQTL_matrix_p6.tsv"
OUTPUT_COVARIATE_FILE="${OUTPUT_DIR}/covariate_matrix_p6.tsv"

# Path to external covariate files (PLACEHOLDER - UPDATE THESE PATHS)
AEI_FILE="./covariates/AEI_Allu_Editing_Index.tsv"
PC_FILE="./covariates/genotype_pcs.tsv"
PEER_FILE="./covariates/PEER_factors.tsv"
CELL_PROP_FILE="./covariates/cell_proportions.tsv"

# --- Setup ---
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Load necessary modules (adjust for your cluster environment)
# module load python/3.9 anaconda/latest # Example modules

# --- Execution ---
echo "Starting Phase 6: Normalization and Covariate Adjustment"
echo "Input Matrix: ${INPUT_FILE}"
echo "Output Normalized Matrix: ${OUTPUT_NORMALIZED_FILE}"
echo "Output Covariate Matrix: ${OUTPUT_COVARIATE_FILE}"
echo "Python Script: ${SCRIPT}"

# Run the Python script
python3 ${SCRIPT} \
    --input_file ${INPUT_FILE} \
    --aei_file ${AEI_FILE} \
    --pc_file ${PC_FILE} \
    --peer_file ${PEER_FILE} \
    --cell_prop_file ${CELL_PROP_FILE} \
    --output_normalized_file ${OUTPUT_NORMALIZED_FILE} \
    --output_covariate_file ${OUTPUT_COVARIATE_FILE}

# --- Success Flag ---
if [ $? -eq 0 ]; then
    echo "Phase 6 completed successfully."
    touch "${OUTPUT_DIR}/phase6_success.flag"
else
    echo "Phase 6 FAILED. Check logs for details."
    exit 1
fi
