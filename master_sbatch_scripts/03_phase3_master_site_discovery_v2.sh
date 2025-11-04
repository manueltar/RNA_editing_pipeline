#!/bin/bash
#
# SLURM Script for Phase 3 (Reworked): Master Site Discovery
# Aggregates all raw calls to create a list of unique, consensus-supported sites.
#

# === SLURM Directives ===
#SBATCH --job-name=P3_Master_Sites
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G                     # High memory for file aggregation
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/p3_master_site_%j.out
#SBATCH --error=logs/p3_master_site_%j.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x

# --- Define Paths ---
DISCOVERY_PYTHON_SCRIPT="./Python_scripts/run_phase3_master_discovery_v2.py"
INPUT_DIR="/scratch/phase2/raw_calls"
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
    --raw_calls_dir ${INPUT_DIR} \
    --output_file ${MASTER_SITE_LIST} \
    --min_edit_level ${MIN_EDIT_LEVEL} # <--- NEW PARAMETER

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 3 Master Site Discovery failed."
    exit 1
fi

echo "Phase 3 complete. Master site list saved to ${MASTER_SITE_LIST}"
