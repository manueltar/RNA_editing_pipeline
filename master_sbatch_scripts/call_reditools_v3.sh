#!/bin/bash
#
# SLURM Script for Phase 2, Step 2: RNA Editing Calling with REDItools3
# V2: Adapted for Array Job processing 250 Chunks and includes Blacklist
#

# === SLURM Directives (Essential for Array Job Execution) ===
#SBATCH --job-name=REDItools_Chunk_%a
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=logs/reditools_call_%A_%a.out
#SBATCH --error=logs/reditools_call_%A_%a.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x
module load samtools/1.18

# --- Configuration (Chunking Logic) ---
CHUNK_MASTER_LIST="/path/to/final/pseudobulks/all_chunk_lists.txt"
INPUT_CHUNK_LIST=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${CHUNK_MASTER_LIST})

# Define the final Python wrapper script path
REDITOOLS_SCRIPT="/path/to/reditools/run_reditools_v5.py" 
REFERENCE_FASTA="/path/to/reference_data/Genome/GRCh38.no_alt.fa"
OUTPUT_DIR="/scratch/phase2/raw_calls"

# --- Blacklist Path (Mandatory Alu/Simple Repeat Filter) ---
BLACKLIST_BED="/path/to/your/reference_data/Annotation/GRCh38_SimpleAlu_Blacklist.bed"

# Create output file path specific to this chunk
CHUNK_ID=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
OUTPUT_FILE="${OUTPUT_DIR}/raw_reditools_output_chunk_${CHUNK_ID}.tsv"

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting REDItools3 calling for Chunk ${SLURM_ARRAY_TASK_ID}..."
echo "Input BAM List: ${INPUT_CHUNK_LIST}"

# --- Execution ---
if [ ! -f "${INPUT_CHUNK_LIST}" ]; then
    echo "CRITICAL ERROR: Chunk list file not found: ${INPUT_CHUNK_LIST}"
    exit 1
fi

# We use the explicit pipeline argument --bed_file to match REDItools3's underlying flag
python3 ${REDITOOLS_SCRIPT} \
    --input_list ${INPUT_CHUNK_LIST} \
    --fasta ${REFERENCE_FASTA} \
    --output ${OUTPUT_FILE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --min_qual 20 \
    --bed_file ${BLACKLIST_BED} # <-- EXPLICITLY MAPPED TO REDITOOLS3 FLAG

if [ $? -ne 0 ]; then
    echo "ERROR: REDItools3 calling failed for Chunk ${CHUNK_ID}."
    exit 1
fi

echo "REDItools3 calling finished. Raw calls saved to ${OUTPUT_FILE}"
 
