#!/bin/bash
# 
# SLURM Script for Phase 2, Step 1: RNA Editing Calling with RED-ML
# V2: Adapted for Array Job processing 250 Chunks of single-cell BAMs
#

# === SLURM Directives ===
#SBATCH --job-name=RED-ML_Chunk_%a
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/redml_call_%A_%a.out
#SBATCH --error=logs/redml_call_%A_%a.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x
module load samtools/1.18 # Required for samtools mpileup/view/etc, if RED-ML uses it.

# --- Configuration ---
# The master list of 250 chunk list paths (from Step 0)
CHUNK_MASTER_LIST="/path/to/final/pseudobulks/all_chunk_lists.txt" 
# This is the path to the text file containing paths of thousands of single-cell BAMs for this chunk
INPUT_CHUNK_LIST=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${CHUNK_MASTER_LIST})

REDML_SCRIPT="/path/to/redml/run_redml_v3.py"
REFERENCE_FASTA="/path/to/reference_data/Genome/GRCh38.no_alt.fa"
OUTPUT_DIR="/scratch/phase2/raw_calls"

# Create output file path specific to this chunk
CHUNK_ID=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
OUTPUT_FILE="${OUTPUT_DIR}/raw_redml_output_chunk_${CHUNK_ID}.tsv"

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting RED-ML calling for Chunk ${SLURM_ARRAY_TASK_ID}..."
echo "Input BAM List: ${INPUT_CHUNK_LIST}"

# --- Execution ---
if [ ! -f "${INPUT_CHUNK_LIST}" ]; then
    echo "CRITICAL ERROR: Chunk list file not found: ${INPUT_CHUNK_LIST}"
    exit 1
fi

# The RED-ML tool must be capable of taking a list of BAMs and processing them as a single group.
# We pass the list of BAMs to RED-ML. (Conceptual command, requires RED-ML support for batch processing)
python3 ${REDML_SCRIPT} \
    --input_list ${INPUT_CHUNK_LIST} \
    --fasta ${REFERENCE_FASTA} \
    --output ${OUTPUT_FILE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --min_qual 20

if [ $? -ne 0 ]; then
    echo "ERROR: RED-ML calling failed for Chunk ${CHUNK_ID}."
    exit 1
fi

echo "RED-ML calling finished. Raw calls saved to ${OUTPUT_FILE}"
