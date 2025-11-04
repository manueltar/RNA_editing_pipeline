#!/bin/bash

# SLURM Directives for the Cleanup Job
# This job should be quick and low-resource, as its only task is file deletion.
#SBATCH --job-name=P7_Cleanup
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00

# --- Get Arguments ---
# Arguments are passed from the master wrapper (000_Wraper_v17.sh)
INDIVIDUAL_ID=$1
INDIVIDUAL_OUTPUT_DIR=$2

# Input validation
if [[ -z "$INDIVIDUAL_ID" || -z "$INDIVIDUAL_OUTPUT_DIR" ]]; then
    echo "ERROR: Missing arguments. Usage: $0 <INDIVIDUAL_ID> <OUTPUT_DIR>"
    exit 1
fi

echo "--- Starting Cleanup for Individual: ${INDIVIDUAL_ID} at $(date) ---"
echo "Targeting output directory: ${INDIVIDUAL_OUTPUT_DIR}"

# --- Intermediate Directories to Remove ---
# These directories are assumed to contain the large, per-chunk/per-caller intermediate 
# files that are no longer needed after Phase 6 is complete.
INTERMEDIATE_DIRS=(
    "P2_REDML_raw"          # Raw output files from REDML for 31 chunks
    "P2_REDItools_raw"      # Raw output files from REDItools for 31 chunks
    "P3_site_chunks"        # Intermediate files from P3 aggregation
    "P4_temp_quant"         # Temporary files created during P4 quantification
)

# --- Perform Cleanup ---
CLEANUP_SUCCESS=true

for dir_name in "${INTERMEDIATE_DIRS[@]}"; do
    TARGET="${INDIVIDUAL_OUTPUT_DIR}/${dir_name}"
    
    if [[ -d "${TARGET}" ]]; then
        echo "Removing intermediate directory: ${TARGET}"
        # Use 'rm -rf' for forced, recursive removal
        rm -rf "${TARGET}"
        
        if [ $? -ne 0 ]; then
            echo "WARNING: Failed to remove directory ${TARGET}. Check permissions."
            CLEANUP_SUCCESS=false
        else
            echo "Successfully removed: ${TARGET}"
        fi
    else
        echo "Intermediate directory not found (Skipping): ${TARGET}"
    fi
done

if $CLEANUP_SUCCESS; then
    echo "--- Cleanup for ${INDIVIDUAL_ID} completed successfully at $(date) ---"
    # Exit 0 so the next individual's job chain can begin (enforced by the 'srun --dependency' in the wrapper)
    exit 0
else
    echo "--- Cleanup for ${INDIVIDUAL_ID} completed with warnings at $(date) ---"
    exit 0 
fi
