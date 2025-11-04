#!/bin/bash
#
# SLURM Script for Phase 3 (Reworked): Individual Site Processing and Normalization
# This script runs once per individual (launched by an outer array job).
# Aggregates raw calls, finds intersection, annotates (using REDIPortal), and normalizes the matrix.
#

# === SLURM Directives ===
#SBATCH --job-name=P3_Proc_Ind
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4         
#SBATCH --mem=16G                 
#SBATCH --time=04:00:00           
#SBATCH --partition=compute
#SBATCH --output=logs/p3_proc_%j.out
#SBATCH --error=logs/p3_proc_%j.err

set -e -o pipefail

# === Arguments from Wrapper (OUTER array job) ===
INDIVIDUAL_ID=$1

# --- Global Reference Paths ---
GLOBAL_REF_DIR="/path/to/your/reference_data/Annotation"
# CRITICAL FIX: Changed from REDIbase to REDIPortal
REDIPortal_BED="${GLOBAL_REF_DIR}/REDIPortal_KnownSites.bed.gz" 
GENE_MAP_TSV="${GLOBAL_REF_DIR}/master_site_gene_map.tsv" 

# --- Define Paths ---
PROCESSING_PYTHON_SCRIPT="./Python_scripts/run_phase3_individual_processing_v4.py"

ROOT_PROJECT_OUTPUT_DIR="/scratch/rna_editing_project/output_data" 

INDIVIDUAL_OUTPUT_DIR="${ROOT_PROJECT_OUTPUT_DIR}/${INDIVIDUAL_ID}/P3_Processed_Matrix"
mkdir -p "${INDIVIDUAL_OUTPUT_DIR}"

FINAL_OUTPUT_FILE="${INDIVIDUAL_OUTPUT_DIR}/${INDIVIDUAL_ID}_processed_matrix.tsv"

# --- Parameters ---
MIN_EDIT_LEVEL=0.1 

# Input validation
if [[ -z "$INDIVIDUAL_ID" ]]; then
    echo "ERROR: Missing INDIVIDUAL_ID argument. Exiting."
    exit 1
fi
if [ ! -f "${REDIPortal_BED}" ] || [ ! -f "${GENE_MAP_TSV}" ]; then
    echo "CRITICAL ERROR: Annotation reference files (REDIPortal or Gene Map) not found. Check Phase 1 setup."
    exit 1
fi


echo "Starting Phase 3 Individual Processing for ${INDIVIDUAL_ID}..."

# --- Execution: Call Python Script ---
python3 ${PROCESSING_PYTHON_SCRIPT} \
    --individual_id "${INDIVIDUAL_ID}" \
    --root_search_dir "${ROOT_PROJECT_OUTPUT_DIR}" \
    --output_file "${FINAL_OUTPUT_FILE}" \
    --min_edit_level ${MIN_EDIT_LEVEL} \
    --rediportal_bed "${REDIPortal_BED}" \
    --gene_map_tsv "${GENE_MAP_TSV}"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 3 Individual Processing failed."
    exit 1
fi

echo "Phase 3 complete. Processed matrix saved to ${FINAL_OUTPUT_FILE}"
