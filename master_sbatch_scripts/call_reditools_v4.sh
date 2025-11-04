#!/bin/bash
#
# SLURM Script for Phase 2: RNA Editing Calling (REDItools3)
# V7: Final Correction based on analyze.py source code. Implements correct QC flags.
#
# === SLURM Directives ===
#SBATCH --job-name=P2_REDItools_Call
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Use 4 cores for mpileup processing
#SBATCH --mem=16G
#SBATCH --time=12:00:00           
#SBATCH --partition=long
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# --- CRITICAL ROBUSTNESS IMPROVEMENT ---
# Exit immediately if a command exits with a non-zero status.
set -e -o pipefail

# === Arguments from Wrapper ===
INDIVIDUAL_ID=$1
INPUT_BAM_ROOT=$2      # e.g., /path/to/input/bams/base_directory/IND_001/cell_type_bams
INDIVIDUAL_OUTPUT_DIR=$3

# --- Global Reference Paths (MUST BE STATIC AND MATCH PHASE 1) ---
GLOBAL_REF_DIR="/path/to/your/reference_data/Annotation"
GENOME_FASTA="/path/to/your/reference_data/Genome/GRCh38.no_alt.fa"

# CRITICAL PATH FIX: Using the combined Master Blacklist file name for exclusion regions
MASTER_BLACKLIST_BED="${GLOBAL_REF_DIR}/GRCh38_SimpleAlu_Master_Blacklist.bed"

# --- TOOL PATHS ---
# REDItools3 is run as a Python module, so the executable path is simplified
# python3 -m reditools analyze [options]
# -----------------

# Input validation
if [[ -z "$INDIVIDUAL_ID" || -z "$INPUT_BAM_ROOT" || -z "$INDIVIDUAL_OUTPUT_DIR" ]]; then
    echo "ERROR: Missing arguments. Exiting array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi
if [ ! -f "${MASTER_BLACKLIST_BED}" ] || [ ! -f "${GENOME_FASTA}" ]; then
    echo "CRITICAL ERROR: Required reference files (FASTA or MASTER_BLACKLIST) not found. Check Phase 1 setup."
    exit 1
fi

# --- Environment Setup ---
module purge
module load python/3.x 
module load samtools/1.18 
# NOTE: Ensure REDItools3 is available in this Python environment.

# --- 1. Dynamic Discovery of Input BAM File (31-Task Logic) ---
echo "Searching for BAM files in: ${INPUT_BAM_ROOT}"
BAM_FILENAMES=($(find "${INPUT_BAM_ROOT}" -maxdepth 1 -type f -name "*.bam" | sort | sed -r 's/.*\/(.*)/\1/'))

BAM_COUNT=${#BAM_FILENAMES[@]}
EXPECTED_COUNT=31

if [ "$BAM_COUNT" -ne "$EXPECTED_COUNT" ]; then
    echo "CRITICAL ERROR: Found ${BAM_COUNT} BAM files. Expected ${EXPECTED_COUNT} (31 cell types)."
    exit 1
fi

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
BAM_FILENAME="${BAM_FILENAMES[${ARRAY_INDEX}]}"
INPUT_BAM_PATH="${INPUT_BAM_ROOT}/${BAM_FILENAME}"
CELL_TYPE_NAME=$(basename "${BAM_FILENAME}" .pseudobulk.bam) 

# --- 2. Define Individual-Specific Output Path ---

OUTPUT_SUBDIR="${INDIVIDUAL_OUTPUT_DIR}/P2_REDItools_raw"
mkdir -p "${OUTPUT_SUBDIR}"

OUTPUT_CALLS_FILE="${OUTPUT_SUBDIR}/${INDIVIDUAL_ID}_${CELL_TYPE_NAME}_reditools_raw.tsv"

echo "--- Starting REDItools3 Call for ${CELL_TYPE_NAME} (${INDIVIDUAL_ID}) ---"
echo "Input BAM: ${INPUT_BAM_PATH}"
echo "Output File: ${OUTPUT_CALLS_FILE}"

# --- 3. Execution: Call REDItools3 (Using correct module command and flags) ---

python3 -m reditools analyze \
    "${INPUT_BAM_PATH}" \
    --reference "${GENOME_FASTA}" \
    --output-file "${OUTPUT_CALLS_FILE}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --min-read-quality 20 \
    --min-base-quality 20 \
    --min-read-depth 10 \
    --min-edits 3 \
    --exclude-regions "${MASTER_BLACKLIST_BED}" 

# --- 4. Final Verification ---

if [ $? -ne 0 ]; then
    echo "ERROR: REDItools3 execution failed. Check log for details."
    exit 1
fi

if [ ! -s "${OUTPUT_CALLS_FILE}" ]; then
    echo "CRITICAL ERROR: REDItools3 output file is missing or empty. This may indicate a runtime failure."
    exit 1
fi

echo "REDItools3 calling successfully completed for ${CELL_TYPE_NAME}."
exit 0
