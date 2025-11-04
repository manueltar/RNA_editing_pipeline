#!/bin/bash
#SBATCH --job-name=Chunk_Prep
#SBATCH --array=1-250%50 # N = total number of chunks (batches)
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --output=logs/chunk_prep_%A_%a.out
 
set -e -o pipefail

# The list of 250 file paths (from Step 0)
CHUNK_MASTER_LIST="/path/to/final/pseudobulks/all_chunk_lists.txt" 
# This variable now holds the path to one of the 250 text files (e.g., /path/to/final/pseudobulks/chunk_list_005)
INPUT_CHUNK_LIST=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${CHUNK_MASTER_LIST})

if [ -f "${INPUT_CHUNK_LIST}" ]; then
    echo "Processing chunk list file: ${INPUT_CHUNK_LIST}"
    # The output of this script is its input, which will be passed to Phase 2 scripts.
    # We simply ensure the path is outputted to standard output for logging, and the file exists.
    # No actual indexing happens here.
else
    echo "CRITICAL ERROR: Chunk list file not found for array ID ${SLURM_ARRAY_TASK_ID}: ${INPUT_CHUNK_LIST}"
    exit 1
fi

# We assume this script's primary role is now to launch the *next* job's array task by passing the chunk list ID.
# Since this script runs as an array task, its success signals that the list exists.

# We must ensure the master wrapper knows the new name of the file to use:
echo "${INPUT_CHUNK_LIST}"
