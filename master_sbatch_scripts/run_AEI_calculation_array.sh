#!/bin/bash
#
# SLURM ARRAY WRAPPER for AEI_calculation: Alu Editing Index Calculation
# Calculates AEI for every individual and every cell type using the Alu-only raw calls.

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_AEI_Master
#SBATCH --output=logs/aei_master_%j.out
#SBATCH --error=logs/aei_master_%j.err
#SBATCH --partition=short        # Master submission script can use 'short'
#SBATCH --mem=1G                 
#SBATCH --time=00:05:00          

# --- Configuration ---
# Directory containing all Phase 2 REDItools output files for AEI
INPUT_DIR="./P2_REDItools_ALU_ONLY" 
# Output directory for individual AEI results
OUTPUT_DIR="./aei_calculation_results"
# BED file containing the coordinates of all Alu elements (Used by reditools index)
ALU_BED_FILE="/path/to/your/reference_data/Annotation/01_ALU_Blacklist.bed" # Synchronized path
# Collation script to run after the array completes
COLLATION_SCRIPT="./Python_scripts/collate_aei_results_p6a.py"
COLLATION_OUTPUT="./phase6_normalized_edQTL/aei_covariate_matrix_p6a.tsv"

# --- Setup ---
mkdir -p ${OUTPUT_DIR}
mkdir -p logs
mkdir -p ./phase6_normalized_edQTL # Ensure the final output dir exists

# Find all gzipped REDItools output files and store their paths
mapfile -t INPUT_FILES < <(find ${INPUT_DIR} -name '*ALU_ONLY.tsv.gz')

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "ERROR: No ALU_ONLY input files found in ${INPUT_DIR}. Please confirm Phase 2 correction was executed."
    exit 1
fi

export ARRAY_MAX_INDEX=$((${#INPUT_FILES[@]} - 1))
echo "Found ${#INPUT_FILES[@]} ALU-ONLY files. Submitting array job 0-${ARRAY_MAX_INDEX}."

# --- Array Job Submission: Part 1 (Calculation) ---
ARRAY_JOB_ID=$(sbatch --parsable --array=0-${ARRAY_MAX_INDEX} \
    --job-name=edQTL_AEI_Calc \
    --output=logs/aei_calculation_%A_%a.out \
    --error=logs/aei_calculation_%A_%a.err \
    --partition=compute \
    --mem=4G \
    --time=00:30:00 \
    --export=ALL,INPUT_FILES,OUTPUT_DIR,ALU_BED_FILE \
    <<EOF
#!/bin/bash
# This is the actual array task script execution

INPUT_FILE=\${INPUT_FILES[\$SLURM_ARRAY_TASK_ID]}
FILENAME=\$(basename "\$INPUT_FILE")
OUTPUT_FILE="${OUTPUT_DIR}/\${FILENAME%.tsv.gz}.aei.tsv"

# REDItools Index Execution: Only counts A->G substitutions within Alu regions
reditools index -B \${ALU_BED_FILE} \${INPUT_FILE} > \${OUTPUT_FILE} 

if [ \$? -ne 0 ]; then
    echo "ERROR: REDItools index failed for \${FILENAME}."
    exit 1
fi
EOF
)

if [ -z "${ARRAY_JOB_ID}" ]; then
    echo "ERROR: Failed to submit AEI Calculation array job. Aborting."
    exit 1
fi
echo "AEI Calculation Array submitted (Job ID: ${ARRAY_JOB_ID})"

# --- Collation Job Submission: Part 2 (Collation) ---
echo "Submitting AEI Collation, dependent on array job completion (afterok:${ARRAY_JOB_ID})..."

sbatch --job-name=edQTL_AEI_Collate \
    --output=logs/aei_collation_%j.out \
    --error=logs/aei_collation_%j.err \
    --partition=short \
    --mem=8G \
    --time=01:00:00 \
    --dependency=afterok:${ARRAY_JOB_ID} \
    --wrap="python3 ${COLLATION_SCRIPT} --input_dir ${OUTPUT_DIR} --output_file ${COLLATION_OUTPUT}"

echo "AEI Collation submitted. AEI_calculation phase is running."
