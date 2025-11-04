#!/bin/bash
#SBATCH --job-name=P2_Cell_Prep
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%A_%a.out

# --- CRITICAL ROBUSTNESS IMPROVEMENT ---
# Exit immediately if a command exits with a non-zero status.
set -e -o pipefail

# --- Arguments from Wrapper ---
# Expected Arguments: $1=INDIVIDUAL_ID, $2=INPUT_BAM_ROOT, $3=INDIVIDUAL_OUTPUT_DIR
INDIVIDUAL_ID=$1
INPUT_BAM_ROOT=$2      # e.g., /path/to/input/bams/base_directory/IND_001/cell_type_bams
INDIVIDUAL_OUTPUT_DIR=$3

# Input validation
if [[ -z "$INDIVIDUAL_ID" || -z "$INPUT_BAM_ROOT" || -z "$INDIVIDUAL_OUTPUT_DIR" ]]; then
    echo "ERROR: Missing arguments. Exiting array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

# --- Environment Setup ---
module purge
module load samtools/1.18

# --- 1. Dynamic Discovery of Input BAM Files ---

# Find all BAM files in the input directory, sorted alphabetically, and store the basenames in a BASH array.
# The `find` command must use the specific input path for the current individual's cell-type BAMs.
# The use of `| sed -r 's/.*\/(.*)/\1/'` extracts only the filename (basename) without the path.

echo "Searching for BAM files in: ${INPUT_BAM_ROOT}"

# CRITICAL FIX: Dynamically populate the list of BAM files
BAM_FILENAMES=($(find "${INPUT_BAM_ROOT}" -maxdepth 1 -type f -name "*.bam" | sort | sed -r 's/.*\/(.*)/\1/'))

# --- 2. Validation of Discovered Files against Array Size ---

# Get the total number of discovered BAMs
BAM_COUNT=${#BAM_FILENAMES[@]}
EXPECTED_COUNT=31 # Hardcoded expectation from the wrapper's ARRAY_SIZE

if [ "$BAM_COUNT" -ne "$EXPECTED_COUNT" ]; then
    echo "CRITICAL ERROR: Found ${BAM_COUNT} BAM files. Expected exactly ${EXPECTED_COUNT} (31 cell types)."
    echo "Please check if the directory ${INPUT_BAM_ROOT} contains the 31 required BAMs."
    exit 1
fi

# --- 3. Determine Current Task's Input File ---

# SLURM array tasks start at 1. Bash arrays start at 0.
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Final check to ensure the array index is valid for the number of files found.
if [ "$ARRAY_INDEX" -ge "$BAM_COUNT" ]; then
    echo "CRITICAL ERROR: SLURM Array ID ${SLURM_ARRAY_TASK_ID} is out of bounds for the ${BAM_COUNT} files found."
    exit 1
fi

BAM_FILENAME="${BAM_FILENAMES[${ARRAY_INDEX}]}"
INPUT_BAM_PATH="${INPUT_BAM_ROOT}/${BAM_FILENAME}"
OUTPUT_BAM_INDEX="${INPUT_BAM_PATH}.bai"

echo "--- Processing Array Task ${SLURM_ARRAY_TASK_ID} for ${INDIVIDUAL_ID} ---"
echo "Target BAM File: ${BAM_FILENAME}"
echo "Full Input Path: ${INPUT_BAM_PATH}"

# --- 4. Core Task: Validate and Index BAM File ---

if [ ! -f "${INPUT_BAM_PATH}" ]; then
    echo "CRITICAL ERROR: Input BAM file not found at ${INPUT_BAM_PATH}. Cannot index. (Should have been caught in file discovery.)"
    exit 1
fi

if [ -f "${OUTPUT_BAM_INDEX}" ]; then
    echo "BAM index (.bai) already exists. Skipping indexing."
else
    echo "Creating BAM index..."
    # Samtools index creates the .bai file directly in the same location
    samtools index "${INPUT_BAM_PATH}"
    
    # Check for index success
    if [ ! -f "${OUTPUT_BAM_INDEX}" ]; then
        echo "CRITICAL ERROR: Samtools indexing failed for ${BAM_FILENAME}."
        exit 1
    fi
fi

# --- Final Check and Output ---
echo "Phase 2 Preparation for ${BAM_FILENAME} complete. BAM is ready for downstream callers."
exit 0
