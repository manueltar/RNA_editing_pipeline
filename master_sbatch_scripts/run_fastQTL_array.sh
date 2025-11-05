#!/bin/bash
#
# SLURM Array Job for Phase 7: FastQTL edQTL Association Testing
#
# This script is designed to be run as an array job, where each task processes
# a specific chunk of the normalized edQTL feature matrix.
#
# CRITICAL: Adjust the array size (1-N) below to match the number of chunks you want.
# Example: 1-100 for 100 chunks.

#SBATCH --job-name=FastQTL_edQTL
#SBATCH --output=logs/fastqtl_%A_%a.out
#SBATCH --error=logs/fastqtl_%A_%a.err
#SBATCH --partition=compute          # Or your preferred partition
#SBATCH --cpus-per-task=4            # Use multiple cores for FastQTL performance
#SBATCH --mem=16G                    # Sufficient memory for one chunk
#SBATCH --time=12:00:00              # Allow plenty of time for association testing
#SBATCH --array=1-100                # <<< ADJUST THIS NUMBER (N=Total Chunks)

# --- Configuration ---
# Set the total number of chunks (must match the --array size)
TOTAL_CHUNKS=100                     # <<< ADJUST THIS NUMBER

# Input/Output Paths
FASTQTL_BIN="/path/to/fastQTL"        # <<< UPDATE THIS PATH
VCF_FILE="./genotype_data/final_genotypes.vcf.gz"  # Your genotype data
NORM_MATRIX="./phase6_normalized_edQTL/normalized_edQTL_matrix_p6.tsv" # Normalized features
COV_MATRIX="./phase6_normalized_edQTL/covariate_matrix_p6.tsv" # Covariates
OUTPUT_DIR="./phase7_edQTL_results"
TEMP_DIR="/scratch/fastqtl_tmp/$SLURM_JOB_ID"

# --- Setup ---
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

TASK_ID=$SLURM_ARRAY_TASK_ID
CHUNK_FILE="${TEMP_DIR}/chunk_${TASK_ID}.bed"
OUTPUT_FILE="${OUTPUT_DIR}/edQTL_results_chunk_${TASK_ID}.tsv"

echo "Starting FastQTL task ${TASK_ID}/${TOTAL_CHUNKS}..."
echo "Feature Matrix: ${NORM_MATRIX}"

# --- 1. Determine Lines to Extract ---
# We need to know the total number of features (lines) in the normalized matrix
# We skip the first line (header) for counting features.
TOTAL_FEATURES=$(cat ${NORM_MATRIX} | tail -n +2 | wc -l)

# Calculate how many features go into each chunk
LINES_PER_CHUNK=$(( (TOTAL_FEATURES + TOTAL_CHUNKS - 1) / TOTAL_CHUNKS ))

# Calculate the starting line (1-based index)
START_LINE=$(( (TASK_ID - 1) * LINES_PER_CHUNK + 2 )) # +2 to skip the main header
END_LINE=$(( TASK_ID * LINES_PER_CHUNK + 1 ))

# Ensure the END_LINE does not exceed the total features + 1 (header)
MAX_END_LINE=$(( TOTAL_FEATURES + 1 ))
if [ $END_LINE -gt $MAX_END_LINE ]; then
    END_LINE=$MAX_END_LINE
fi

NUM_LINES_TO_EXTRACT=$(( END_LINE - START_LINE + 1 ))

echo "Total Features: ${TOTAL_FEATURES}"
echo "Lines per chunk: ${LINES_PER_CHUNK}"
echo "Extracting lines ${START_LINE} through ${END_LINE} (${NUM_LINES_TO_EXTRACT} lines)..."

# --- 2. Extract and Format the Chunk ---
# FastQTL requires the header in the first line, followed by the features.
# It also strictly expects the first 4 columns to be Chromosome, Start, End, and FeatureID.
# Since your Phase 6 output is TSV, we'll extract the full lines, including the header for the first chunk.

# Extract the header (first line)
head -n 1 ${NORM_MATRIX} > ${CHUNK_FILE}

# Extract the feature rows for this chunk
tail -n +${START_LINE} ${NORM_MATRIX} | head -n ${NUM_LINES_TO_EXTRACT} >> ${CHUNK_FILE}

# --- 3. Run FastQTL ---
# Assuming FastQTL takes normalized feature matrix (-bed), genotype VCF (-vcf), and covariates (-cov)
echo "Running FastQTL..."
${FASTQTL_BIN} \
    --vcf ${VCF_FILE} \
    --bed ${CHUNK_FILE} \
    --cov ${COV_MATRIX} \
    --out ${OUTPUT_FILE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --chunk ${TASK_ID} ${TOTAL_CHUNKS} # FastQTL chunking utility

# --- 4. Cleanup and Validation ---
if [ $? -eq 0 ]; then
    echo "FastQTL job ${TASK_ID} completed successfully. Results saved to ${OUTPUT_FILE}"
    rm ${CHUNK_FILE} # Clean up the temporary chunk file
else
    echo "ERROR: FastQTL job ${TASK_ID} failed."
    exit 1
fi
