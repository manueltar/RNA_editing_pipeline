#!/bin/bash
#
# SLURM Script for Phase 6: Trait Normalization (Inverse Normal Transformation)
# Applies filtering and INT to the gene-level AvgER matrix, preparing the final trait matrix for edQTL mapping.
#

# === SLURM Directives ===
#SBATCH --job-name=P6_Normalization
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/p6_normalization_%j.out
#SBATCH --error=logs/p6_normalization_%j.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x
# Note: R/SciPy is needed for the Inverse Normal Transformation

# --- Define Paths ---
NORMALIZATION_PYTHON_SCRIPT="./Python_scripts/run_phase6_normalization.py"

# --- Input/Output Files ---
# Input from Phase 5 (Gene Summary)
INPUT_MATRIX="/scratch/phase5/final_summary/gene_editing_summary_p5.tsv"

# Final Output (Normalized Trait Matrix)
OUTPUT_DIR="/scratch/phase6/normalized_traits"
NORMALIZED_TRAIT_MATRIX="${OUTPUT_DIR}/normalized_avg_editing_traits_p6.tsv"

# --- Parameters ---
# Minimum number of samples (cell types) in which a gene must have quantifiable editing 
# (i.e., less than 50% NA values for the AvgER to be reliable, assuming 31 cell types/samples)
MIN_QUANTIFIABLE_SAMPLES=15 

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting Phase 6: Trait Normalization and Filtering..."

if [ ! -f "${INPUT_MATRIX}" ]; then
    echo "ERROR: Input matrix from Phase 5 not found: ${INPUT_MATRIX}"
    exit 1
fi

# --- Execution: Call Python Script ---
python3 ${NORMALIZATION_PYTHON_SCRIPT} \
    --input_matrix ${INPUT_MATRIX} \
    --output_file ${NORMALIZED_TRAIT_MATRIX} \
    --min_quantifiable_samples ${MIN_QUANTIFIABLE_SAMPLES}

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 6 Normalization failed."
    exit 1
fi

echo "Phase 6 complete. Normalized trait matrix saved to ${NORMALIZED_TRAIT_MATRIX}"
