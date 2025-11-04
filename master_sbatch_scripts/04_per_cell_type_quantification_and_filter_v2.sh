#!/bin/bash
#
# SLURM Script for Phase 4 (Reworked) v2: Direct input of Phase 3 Matrix
# Queries the 31 BAMs using the full Phase 3 matrix, applies final filters, and annotates.
#

# === SLURM Directives ===
#SBATCH --job-name=P4_Quant_Filter_v2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16          
#SBATCH --mem=48G
#SBATCH --time=12:00:00             
#SBATCH --partition=compute
#SBATCH --output=logs/p4_quant_v2_%j.out
#SBATCH --error=logs/p4_quant_v2_%j.err

set -e -o pipefail

# === Arguments from Wrapper (OUTER array job) ===
INDIVIDUAL_ID=$1 # <--- Assume INDIVIDUAL_ID is passed as the first argument

# --- Define Paths ---
QUANT_PYTHON_SCRIPT="./Python_scripts/run_phase4_quantification_v2.py"

# --- Input/Output Directories ---
ROOT_PROJECT_OUTPUT_DIR="/scratch/rna_editing_project/output_data"
PHASE3_INPUT_DIR="${ROOT_PROJECT_OUTPUT_DIR}/${INDIVIDUAL_ID}/P3_Processed_Matrix"
OUTPUT_DIR="/scratch/rna_editing_project/phase4/quantification_matrix"

# CRITICAL FIX: Direct path to the Phase 3 output matrix
MASTER_SITE_MATRIX="${PHASE3_INPUT_DIR}/${INDIVIDUAL_ID}_annotated_raw_matrix.tsv"

# Reference Data
GERMLINE_VCF="/path/to/individual/germline.vcf.gz" 
GTF_ANNOTATION="/path/to/annotation/Homo_sapiens.GRCh38.109.gtf.gz" 

# Input BAMs directory (31 cell-type BAMs for this individual)
INPUT_BAMS_DIR="/path/to/31_cell_type_bams_for_${INDIVIDUAL_ID}"
BAM_FILE_PATTERN="*sorted.bam" 

# Final Output Matrix
QUANTIFICATION_MATRIX="${OUTPUT_DIR}/${INDIVIDUAL_ID}_final_editing_matrix_p4.tsv"

# --- Parameters (QC and Filtering) ---
SPLICE_SITE_THRESHOLD=4
MIN_READ_COVERAGE=10 
THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting Phase 4 (v2): Quantification and Filtering for ${INDIVIDUAL_ID}..."

if [ ! -f "${MASTER_SITE_MATRIX}" ]; then
    echo "ERROR: Phase 3 Annotated Matrix not found: ${MASTER_SITE_MATRIX}. Cannot proceed."
    exit 1
fi

# --- Execution: Call Python Script ---
python3 ${QUANT_PYTHON_SCRIPT} \
    --phase3_matrix "${MASTER_SITE_MATRIX}" \
    --output_file "${QUANTIFICATION_MATRIX}" \
    --input_bams_dir "${INPUT_BAMS_DIR}" \
    --bam_pattern "${BAM_FILE_PATTERN}" \
    --germline_vcf "${GERMLINE_VCF}" \
    --individual_id "${INDIVIDUAL_ID}" \
    --gtf_annotation "${GTF_ANNOTATION}" \
    --splice_site_threshold ${SPLICE_SITE_THRESHOLD} \
    --min_read_coverage ${MIN_READ_COVERAGE} \
    --threads ${THREADS}

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 4 Quantification failed."
    exit 1
fi

echo "Phase 4 (v2) complete. Quantification matrix saved to ${QUANTIFICATION_MATRIX}"
