#!/bin/bash
#
# SLURM Script for Phase 2: RNA Editing Calling (RED-ML)
# V6: Final Correction for Perl Tool Execution and Correct Output Filename
#
# === SLURM Directives ===
#SBATCH --job-name=P2_REDML_Call
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
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

# Paths to the deconstructed blacklists
SIMPLE_REPEATS_BED="${GLOBAL_REF_DIR}/01_REPEATS_Blacklist.bed" # For --simpleRepeat
ALU_BED="${GLOBAL_REF_DIR}/01_ALU_Blacklist.bed"                 # For --alu

# --- TOOL PATHS ---
# !!! CRITICAL: DEPLOYABLE LOGIC REQUIRES SETTING THIS PATH !!!
REDML_EXECUTABLE="/path/to/RED-ML/red_ML.pl"
# -----------------

# Input validation
if [[ -z "$INDIVIDUAL_ID" || -z "$INPUT_BAM_ROOT" || -z "$INDIVIDUAL_OUTPUT_DIR" ]]; then
    echo "ERROR: Missing arguments. Exiting array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi
if [ ! -f "${SIMPLE_REPEATS_BED}" ] || [ ! -f "${ALU_BED}" ] || [ ! -f "${GENOME_FASTA}" ]; then
    echo "CRITICAL ERROR: Required reference files for RED-ML (REPEATS/ALU/FASTA) not found. Check Phase 1 setup."
    exit 1
fi

# --- Environment Setup ---
module purge
module load perl/5.x 
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

# --- 2. Define Individual-Specific Output Path ---

OUTPUT_SUBDIR="${INDIVIDUAL_OUTPUT_DIR}/P2_REDML_raw"
mkdir -p "${OUTPUT_SUBDIR}"

OUTPUT_CALLS_DIR="${OUTPUT_SUBDIR}/${INDIVIDUAL_ID}_${CELL_TYPE_NAME}_redml_raw_outdir"
mkdir -p "${OUTPUT_CALLS_DIR}"

echo "--- Starting RED-ML Call for ${CELL_TYPE_NAME} (${INDIVIDUAL_ID}) ---"

# --- 3. Execution: Call RED-ML (Perl command) ---

perl ${REDML_EXECUTABLE} \
    --rnabam "${INPUT_BAM_PATH}" \
    --reference "${GENOME_FASTA}" \
    --simpleRepeat "${SIMPLE_REPEATS_BED}" \
    --alu "${ALU_BED}" \
    --outdir "${OUTPUT_CALLS_DIR}" \
    --p 0.5 \
    --threads ${SLURM_CPUS_PER_TASK} 

# --- 4. Final Verification and Renaming ---
# CRITICAL FIX: Using the correct output file name: RNA_editing.sites.txt
FINAL_RAW_CALLS="${OUTPUT_CALLS_DIR}/RNA_editing.sites.txt" 
FINAL_RENAMED_FILE="${OUTPUT_SUBDIR}/${INDIVIDUAL_ID}_${CELL_TYPE_NAME}_redml_raw.tsv"

if [ $? -ne 0 ]; then
    echo "ERROR: RED-ML execution failed. Check log for details."
    exit 1
fi

if [ ! -s "${FINAL_RAW_CALLS}" ]; then
    echo "CRITICAL ERROR: RED-ML output file (${FINAL_RAW_CALLS}) is missing or empty. This may indicate a runtime failure."
    # If the file is missing, we check for the diagnostic file to ensure the tool at least ran
    if [ ! -d "${OUTPUT_CALLS_DIR}" ]; then
        echo "FATAL: Output directory ${OUTPUT_CALLS_DIR} not found. Tool failed early."
    fi
    exit 1
fi

# Rename the output file to the standardized name expected by Phase 3
mv "${FINAL_RAW_CALLS}" "${FINAL_RENAMED_FILE}"

echo "RED-ML calling successfully completed for ${CELL_TYPE_NAME}."
echo "Raw calls saved to ${FINAL_RENAMED_FILE}"
exit 0
