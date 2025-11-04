#!/bin/bash
#
# SLURM Script for Phase 4 (Reworked): Per-Cell Quantification & Final Filter
# Queries the 31 BAMs using the master site list, applies final filters, and annotates.
#

# === SLURM Directives ===
#SBATCH --job-name=P4_Quant_Filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16             
#SBATCH --mem=48G
#SBATCH --time=12:00:00                
#SBATCH --partition=compute
#SBATCH --output=logs/p4_quant_%j.out
#SBATCH --error=logs/p4_quant_%j.err

set -e -o pipefail

# === Setup ===
module purge
module load python/3.x
module load pysam # Ensure pysam environment is available

# --- Define Paths ---
QUANT_PYTHON_SCRIPT="./Python_scripts/run_phase4_quantification.py"

# --- Input/Output Files ---
MASTER_LIST_DIR="/scratch/phase3/master_site_list"
OUTPUT_DIR="/scratch/phase4/quantification_matrix"

# Input from Phase 3
MASTER_SITE_LIST="${MASTER_LIST_DIR}/consensus_master_sites_p3.tsv"

# Reference Data
GERMLINE_VCF="/path/to/individual/germline.vcf.gz"
GTF_ANNOTATION="/path/to/annotation/genes.gtf" 

# Input BAMs directory (The 31 files)
INPUT_BAMS_DIR="/path/to/31_individual_bams"
BAM_FILE_PATTERN="scRNA_file_*.bam" 

# Final Output Matrix
QUANTIFICATION_MATRIX="${OUTPUT_DIR}/final_editing_matrix_p4.tsv"

# --- Parameters (NEW/MODIFIED) ---
INDIVIDUAL_ID="IndA"                  
SPLICE_SITE_THRESHOLD=4
MIN_READ_COVERAGE=10  # <--- New Parameter Defined Here

mkdir -p logs
mkdir -p ${OUTPUT_DIR}

echo "Starting Phase 4: Per-Cell Quantification and Filtering..."

if [ ! -f "${MASTER_SITE_LIST}" ]; then
    echo "ERROR: Master Site List from Phase 3 not found: ${MASTER_SITE_LIST}"
    exit 1
fi

# --- Execution: Call Python Script ---
python3 ${QUANT_PYTHON_SCRIPT} \
    --master_sites ${MASTER_SITE_LIST} \
    --output_file ${QUANTIFICATION_MATRIX} \
    --input_bams_dir ${INPUT_BAMS_DIR} \
    --bam_pattern "${BAM_FILE_PATTERN}" \
    --germline_vcf ${GERMLINE_VCF} \
    --individual_id ${INDIVIDUAL_ID} \
    --gtf_annotation ${GTF_ANNOTATION} \
    --splice_site_threshold ${SPLICE_SITE_THRESHOLD} \
    --min_read_coverage ${MIN_READ_COVERAGE} \
    --threads ${SLURM_CPUS_PER_TASK}

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 4 Quantification failed."
    exit 1
fi

echo "Phase 4 complete. Quantification matrix saved to ${QUANTIFICATION_MATRIX}"
