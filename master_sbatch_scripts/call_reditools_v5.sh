#!/bin/bash
#
# SLURM Script for Phase 2: RNA Editing Calling (REDItools3)
# V9: Adds a second, parallel run of REDItools to generate Alu-only raw calls 
#     required for subsequent Alu Editing Index (AEI) calculation.
#     Paths synchronized with Phase 1 output names.
#
# === SLURM Directives ===
#SBATCH --job-name=P2_REDItools_Call
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Use 4 cores for mpileup processing
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
INPUT_BAM_ROOT=$2      # e.g., /path/to/input/bams/base_directory/IND_001/cell_type_bams
INDIVIDUAL_OUTPUT_DIR=$3

# --- Global Reference Paths (MUST BE STATIC AND MATCH PHASE 1) ---
GLOBAL_REF_DIR="/path/to/your/reference_data/Annotation"
GENOME_FASTA="/path/to/your/reference_data/Genome/GRCh38.no_alt.fa"

# CRITICAL PATH FIX: Master Blacklist file name for exclusion regions (from Phase 1)
MASTER_BLACKLIST_BED="${GLOBAL_REF_DIR}/GRCh38_SimpleAlu_Master_Blacklist.bed"

# NEW CRITICAL PATH: Alu-Only BED file (Subset of the master blacklist)
# Filename corrected to match the Phase 1 script output: 01_ALU_Blacklist.bed
ALU_ONLY_BED="${GLOBAL_REF_DIR}/01_ALU_Blacklist.bed" 

# --- TOOL PATHS ---
# python3 -m reditools analyze [options]
# -----------------

# Input validation
if [[ -z "$INDIVIDUAL_ID" || -z "$INPUT_BAM_ROOT" || -z "$INDIVIDUAL_OUTPUT_DIR" ]]; then
    echo "ERROR: Missing arguments. Exiting array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi
if [ ! -f "${MASTER_BLACKLIST_BED}" ] || [ ! -f "${GENOME_FASTA}" ] || [ ! -f "${ALU_ONLY_BED}" ]; then
    echo "CRITICAL ERROR: Required reference files (FASTA, MASTER_BLACKLIST, or ALU_ONLY_BED) not found. Check Phase 1 setup."
    exit 1
fi

# --- Environment Setup ---
module purge
module load python/3.x 
module load samtools/1.18 

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

# --- 2. Define Individual-Specific Output Paths and Create Directories ---

# Path for standard, non-Alu-filtered output (Input for Phase 3+)
OUTPUT_SUBDIR="${INDIVIDUAL_OUTPUT_DIR}/P2_REDItools_raw"
mkdir -p "${OUTPUT_SUBDIR}"
OUTPUT_CALLS_FILE="${OUTPUT_SUBDIR}/${INDIVIDUAL_ID}_${CELL_TYPE_NAME}_reditools_raw.tsv"

# Path for ALU-ONLY output (Input for AEI_calculation)
AEI_OUTPUT_SUBDIR="${INDIVIDUAL_OUTPUT_DIR}/P2_REDItools_ALU_ONLY"
mkdir -p "${AEI_OUTPUT_SUBDIR}"
AEI_CALLS_FILE="${AEI_OUTPUT_SUBDIR}/${INDIVIDUAL_ID}_${CELL_TYPE_NAME}_ALU_ONLY.tsv"


echo "--- Starting Dual REDItools3 Call for ${CELL_TYPE_NAME} (${INDIVIDUAL_ID}) ---"
echo "Input BAM: ${INPUT_BAM_PATH}"

# --- 3A. Execution: Call REDItools3 for Main edQTL Sites (Excluding Alu) ---
echo "Running REDItools (MAIN edQTL: Excluding Master Blacklist)... Output: ${OUTPUT_CALLS_FILE}"

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


# --- 3B. Execution: Call REDItools3 for AEI Calculation (Alu ONLY) ---
echo "Running REDItools (AEI: Alu ONLY)... Output: ${AEI_CALLS_FILE}"

python3 -m reditools analyze \
    "${INPUT_BAM_PATH}" \
    --reference "${GENOME_FASTA}" \
    --output-file "${AEI_CALLS_FILE}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --min-read-quality 20 \
    --min-base-quality 20 \
    --min-read-depth 10 \
    --min-edits 3 \
    --region "${ALU_ONLY_BED}" # <--- CRITICAL: --region only includes sites in the BED file


# --- 4. Final Verification and Gzip (Updated) ---

if [ $? -ne 0 ]; then
    echo "ERROR: One of the REDItools executions failed. Check log for details."
    exit 1
fi

if [ ! -s "${OUTPUT_CALLS_FILE}" ] || [ ! -s "${AEI_CALLS_FILE}" ]; then
    echo "CRITICAL ERROR: One or both REDItools output files are missing or empty."
    exit 1
fi

echo "Gzipping main output file..."
gzip -f "${OUTPUT_CALLS_FILE}"

echo "Gzipping AEI output file..."
gzip -f "${AEI_CALLS_FILE}"


echo "REDItools3 calling successfully completed for ${CELL_TYPE_NAME}. Both main and AEI files generated."
exit 0
