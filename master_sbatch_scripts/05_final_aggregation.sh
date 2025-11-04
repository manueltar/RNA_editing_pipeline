#!/bin/bash
#
# SLURM Script for Phase 5: Final Aggregation
# Collapses the Phase 4 site-by-cell quantification matrix into a gene-by-cell summary matrix.
#

# === SLURM Directives ===
#SBATCH --job-name=P5_Aggregation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/p5_aggregation_%j.out
#SBATCH --error=logs/p5_aggregation_%j.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x

# --- Define Paths ---
# Ensure this path is correct for your environment
AGGREGATION_PYTHON_SCRIPT="./Python_scripts/run_phase5_aggregation.py"

# --- Input/Output Files ---
# Input from Phase 4
QUANTIFICATION_MATRIX="/scratch/phase4/quantification_matrix/final_editing_matrix_p4.tsv"

# Final Output Summary
OUTPUT_DIR="/scratch/phase5/final_summary"
GENE_SUMMARY_MATRIX="${OUTPUT_DIR}/gene_editing_summary_p5.tsv"

# --- Parameters ---
GROUPING_COLUMN="gene_name" 

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting Phase 5: Final Aggregation..."

if [ ! -f "${QUANTIFICATION_MATRIX}" ]; then
    echo "ERROR: Input matrix from Phase 4 not found: ${QUANTIFICATION_MATRIX}"
    exit 1
fi

# --- Execution: Call Python Script ---
python3 ${AGGREGATION_PYTHON_SCRIPT} \
    --input_matrix ${QUANTIFICATION_MATRIX} \
    --output_file ${GENE_SUMMARY_MATRIX} \
    --group_by ${GROUPING_COLUMN}

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 5 Aggregation failed."
    exit 1
fi

echo "Phase 5 complete. Gene summary matrix saved to ${GENE_SUMMARY_MATRIX}"
