#!/bin/bash
#SBATCH --job-name=phase5_collation
#SBATCH --mem=128G         # IMPORTANT: Increased memory (128GB) is necessary for loading 6000 files into a single matrix. Adjust this if your cluster needs more or less.
#SBATCH --time=06:00:00    # Increased time limit for heavy I/O and processing
#SBATCH --cpus-per-task=8  # Use multiple cores for faster processing
#SBATCH --output=logs/phase5_collation_%j.out
#SBATCH --error=logs/phase5_collation_%j.err
#SBATCH --partition=compute # Use an appropriate partition if your cluster requires it

# --- Configuration ---
# Directory where all 6000 individual Phase 4 files are located
INPUT_DIR="./phase4_output_matrices"
# Directory where the single final Phase 5 matrix will be saved
OUTPUT_DIR="./phase5_eQTL_features"
# The Python script to execute (which handles the collation and selection)
SCRIPT="./collate_and_select_phase5.py"

# Define the final output file name
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/eQTL_features_all_genes_all_celltypes_p5.tsv"

# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p logs

# --- Execution ---
echo "--- Starting Phase 5: Population Collation and Feature Selection ---"
echo "Reading files from: ${INPUT_DIR}"
echo "Output will be saved to: ${FINAL_OUTPUT_FILE}"

# Activate your Python environment (e.g., conda or virtual environment) if needed
# module load python/3.9 # Example for modules
# source /path/to/your/env/bin/activate 

python3 ${SCRIPT} \
    --input_dir ${INPUT_DIR} \
    --file_pattern "*_final_editing_matrix_p4.tsv" \
    --output_file ${FINAL_OUTPUT_FILE}

if [ $? -eq 0 ]; then
    echo "SUCCESS: Phase 5 Collation completed."
    echo "The final eQTL feature matrix is: ${FINAL_OUTPUT_FILE}"
    # Creates a flag file for the master pipeline script to check before proceeding to Phase 6
    touch ${OUTPUT_DIR}/phase5_success.flag
else
    echo "ERROR: Phase 5 Collation failed. Check logs/phase5_collation_${SLURM_JOB_ID}.err"
fi
