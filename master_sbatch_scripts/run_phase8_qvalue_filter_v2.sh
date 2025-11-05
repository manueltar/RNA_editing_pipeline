#!/bin/bash
#
# SLURM wrapper for Phase 8: FDR Correction and Lead SNP Identification
# This script applies the Benjamini/Hochberg procedure to control FDR
# for both edQTL and AEI-QTL mapping results.

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P8_Qvalue
#SBATCH --output=logs/phase8_qvalue_%j.out
#SBATCH --error=logs/phase8_qvalue_%j.err
#SBATCH --partition=bigmem        
#SBATCH --mem=32G                
#SBATCH --time=08:00:00          
#SBATCH --cpus-per-task=4         
#SBATCH --nodes=1

# --- Configuration ---
SCRIPT="./Python_scripts/qvalue_filter_phase8.py"
OUTPUT_DIR="./phase8_final_results"
mkdir -p ${OUTPUT_DIR}

# 1. edQTL Results (Phase 7)
INPUT_EDQTL_DIR="./phase7_edQTL_results"
# NOTE: The AEI-QTL is in a single file, the edQTL needs a list of all chromosome files.
# We assume the main edQTL mapping results are compressed and follow a simple naming pattern.
EDQTL_FILES=$(find ${INPUT_EDQTL_DIR} -name '*fastqtl_raw.txt.gz' -type f)

# 2. AEI-QTL Results (Phase 7b)
INPUT_AEIQTL_DIR="./phase7b_aeiqtl_results"
AEIQTL_FILES=$(find ${INPUT_AEIQTL_DIR} -name 'aei_qtl_mapping_results_p7b.tsv.gz' -type f)

# Check for required input files
if [ -z "$EDQTL_FILES" ] || [ -z "$AEIQTL_FILES" ]; then
    echo "CRITICAL ERROR: FastQTL results from Phase 7 or 7b not found. Check directories."
    exit 1
fi

# --- Execution ---
echo "Starting Phase 8: FDR Correction and Lead SNP Identification"

# Use the list of files found by 'find' as the arguments for nargs='+' in the Python script
python3 ${SCRIPT} \
    --input_edqtl_results ${EDQTL_FILES} \
    --input_aeiqtl_results ${AEIQTL_FILES} \
    --output_dir ${OUTPUT_DIR}

# --- Success Flag ---
if [ $? -eq 0 ]; then
    echo "Phase 8 completed successfully. Final significant results saved to ${OUTPUT_DIR}."
    touch "${OUTPUT_DIR}/phase8_success.flag"
else
    echo "Phase 8 FAILED. Check logs for details."
    exit 1
fi
