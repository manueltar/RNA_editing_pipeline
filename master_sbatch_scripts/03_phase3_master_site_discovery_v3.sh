#!/bin/bash
#
# SLURM Script for Phase 3: Master Site Discovery
# Aggregates all raw calls to create a list of unique, consensus-supported sites.
#

# === SLURM Directives ===
#SBATCH --job-name=P3_Master_Sites
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4         # Increased for potential parallel I/O in Python
#SBATCH --mem=64G                 # High memory for massive file aggregation
#SBATCH --time=08:00:00           # Increased time limit for safety
#SBATCH --partition=compute
#SBATCH --output=logs/p3_master_site_%j.out
#SBATCH --error=logs/p3_master_site_%j.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x

# --- Define Paths ---
DISCOVERY_PYTHON_SCRIPT="./Python_scripts/run_phase3_master_discovery_v3.py"

# CRITICAL FIX: The Python script MUST recursively search the root of all individual outputs
# We assume a common project root directory where all IND_XXXX subdirectories reside.
ROOT_PROJECT_OUTPUT_DIR="/scratch/rna_editing_project/output_data" 

OUTPUT_DIR="/scratch/phase3/master_site_list"

# Final Output (Only coordinates and basic attributes)
MASTER_SITE_LIST="${OUTPUT_DIR}/consensus_master_sites_p3.tsv"

# --- Parameters (NEW) ---
MIN_EDIT_LEVEL=0.1 # Explicit minimum edit level for high-confidence sites

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting Phase 3: Master Site Discovery..."

# --- Execution: Call Python Script ---
python3 ${DISCOVERY_PYTHON_SCRIPT} \
    --root_search_dir ${ROOT_PROJECT_OUTPUT_DIR} \
    --output_file ${MASTER_SITE_LIST} \
    --min_edit_level ${MIN_EDIT_LEVEL} \
    --threads ${SLURM_CPUS_PER_TASK} # Pass threads for file I/O optimization

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 3 Master Site Discovery failed."
    exit 1
fi

echo "Phase 3 complete. Master site list saved to ${MASTER_SITE_LIST}"
