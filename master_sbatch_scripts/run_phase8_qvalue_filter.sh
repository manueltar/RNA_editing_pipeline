#!/bin/bash
#
# SLURM wrapper for Phase 8: Multiple Testing Correction and Filtering
# This script applies FDR correction (Q-value) to the Phase 7 FastQTL results 
# and extracts the final list of significant lead edQTLs.

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P8_Filter
#SBATCH --output=logs/phase8_qvalue_%j.out
#SBATCH --error=logs/phase8_qvalue_%j.err
#SBATCH --partition=short        # Should be a fast, standard partition
#SBATCH --mem=8G                 # Sufficient memory for pandas/numpy processing
#SBATCH --time=02:00:00          # Should complete quickly

# --- Configuration ---
SCRIPT="./Python_scripts/process_fastqtl_results_p8.py"

# Input file path (Output from Phase 7, must be unzipped if FastQTL output was compressed)
INPUT_DIR="./phase7_edqtl_results"
INPUT_FILE="${INPUT_DIR}/edQTL_mapping_results_p7.tsv" 

# Output directory
OUTPUT_DIR="./phase8_final_edqtl_signals"
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/significant_lead_edqtl_p8.tsv"

# Filtering Parameter
FDR_THRESHOLD=0.05

# --- Setup ---
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Ensure input file is uncompressed if necessary (FastQTL output is often compressed)
if [ -f "${INPUT_FILE}.gz" ] && [ ! -f "${INPUT_FILE}" ]; then
    echo "Decompressing FastQTL results..."
    gunzip "${INPUT_FILE}.gz"
    if [ $? -ne 0 ]; then
        echo "Decompression FAILED."
        exit 1
    fi
fi

# Load necessary modules
# module load python/3.9 anaconda/latest # Example modules

# --- Execution ---
echo "Starting Phase 8: Multiple Testing Correction and Filtering"
echo "FDR Threshold: ${FDR_THRESHOLD}"
echo "Input File: ${INPUT_FILE}"

python3 ${SCRIPT} \
    --input_file ${INPUT_FILE} \
    --output_file ${FINAL_OUTPUT_FILE} \
    --fdr_threshold ${FDR_THRESHOLD}

# --- Success Flag ---
if [ $? -eq 0 ]; then
    echo "Phase 8 completed successfully."
    touch "${OUTPUT_DIR}/phase8_success.flag"
else
    echo "Phase 8 FAILED. Check logs for details."
    exit 1
fi
