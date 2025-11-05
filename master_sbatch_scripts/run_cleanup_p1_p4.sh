#!/bin/bash
#
# SLURM wrapper for Cleanup Phase (P1-P4)
# Erases intermediate files to free up disk space before starting P5/Statistical Analysis.

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P1_P4_Cleanup
#SBATCH --output=logs/cleanup_%j.out
#SBATCH --error=logs/cleanup_%j.err
#SBATCH --partition=short
#SBATCH --mem=1G
#SBATCH --time=00:15:00 

# --- Configuration ---
# Directories containing large intermediate files from P1-P4
INTERMEDIATE_DIRS=(
    "./P1_raw_REDItools_output"
    "./P3_REDML_unfiltered_calls"
    "./logs/reditools_raw_calls"
    "./P4_intermediates"
)

# Files/Patterns to Keep (These are the inputs for P5 and AEI-Calc)
# 1. The final, corrected, ALU-ONLY raw call files.
ALU_ONLY_DIR="./P2_REDItools_ALU_ONLY"

# --- Execution ---
echo "--- Starting Cleanup of Intermediate Files (P1-P4) ---"

# 1. Delete intermediate directories
for DIR in "${INTERMEDIATE_DIRS[@]}"; do
    if [ -d "$DIR" ]; then
        echo "Deleting intermediate directory: $DIR"
        rm -rf "$DIR"
    else
        echo "Directory not found (Skipping): $DIR"
    fi
done

# 2. Delete specific unwanted files (RED-ML artifacts)
echo "Deleting specific unwanted RED-ML files..."
find . -type f -name 'variation.sites.feature.txt' -delete -print
find . -type f -name 'mut.txt.gz' -delete -print

# 3. Final check: Ensure the critical input directory is preserved
if [ -d "$ALU_ONLY_DIR" ]; then
    echo "SUCCESS: Preserved critical input directory for P5/AEI-Calc: $ALU_ONLY_DIR"
else
    echo "CRITICAL ERROR: ALU_ONLY input directory is missing. ABORTING CLEANUP."
    exit 1
fi

echo "--- Cleanup Complete. Space freed. ---"
