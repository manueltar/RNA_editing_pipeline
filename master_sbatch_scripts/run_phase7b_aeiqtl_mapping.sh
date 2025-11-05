#!/bin/bash
#
# SLURM wrapper for Phase 7b: AEI-QTL Association Mapping (FastQTL)
# Maps genetic variants to the Alu Editing Index (AEI).

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P7b_AEIQTL
#SBATCH --output=logs/phase7b_aeiqtl_%j.out
#SBATCH --error=logs/phase7b_aeiqtl_%j.err
#SBATCH --partition=bigmem        
#SBATCH --mem=64G                
#SBATCH --time=48:00:00          
#SBATCH --cpus-per-task=16        
#SBATCH --nodes=1

# --- Configuration: Paths ---
PREP_SCRIPT="./Python_scripts/prepare_aeiqtl_matrices.py"
FASTQTL_BIN="/path/to/FastQTL/executable"

# Input files (Outputs from Phase 6 and AEI_calculation)
INPUT_AEI_COVS="./phase6_normalized_edQTL/aei_covariate_matrix_p6a.tsv"
INPUT_COVS_DIR="./input_covariates"
INPUT_GENOTYPE_PCS="${INPUT_COVS_DIR}/genotype_pcs.tsv"
INPUT_PEER_FACTORS="${INPUT_COVS_DIR}/peer_factors_k60.tsv"

# Genotype and Output Paths
GENOTYPE_FILE="./genotype_data/genotype_all_individuals.vcf.gz"
OUTPUT_DIR="./phase7b_aeiqtl_results"
mkdir -p ${OUTPUT_DIR}

# AEI-QTL Specific Files
AEIQTL_PHENOTYPE="${OUTPUT_DIR}/aeiqtl_phenotype_p7b.tsv"
AEIQTL_COVARIATES="${OUTPUT_DIR}/aeiqtl_covariate_matrix_p7b.tsv"
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/aei_qtl_mapping_results_p7b.tsv.gz"
LOG_FILE="${OUTPUT_DIR}/aeiqtl_run_log.txt"

# Mapping Parameters (Same as P7)
CIS_WINDOW_SIZE=1000000
PERMUTATIONS=1000
# --- Setup ---
echo "--- Starting Phase 7b: AEI-QTL Mapping Preparation ---"
date | tee ${LOG_FILE}

# 1. Preparation (Run Python script)
python3 ${PREP_SCRIPT} \
    --input_aei_covariates ${INPUT_AEI_COVS} \
    --input_genotype_pcs ${INPUT_GENOTYPE_PCS} \
    --input_peer_factors ${INPUT_PEER_FACTORS} \
    --output_phenotype ${AEIQTL_PHENOTYPE} \
    --output_covariates ${AEIQTL_COVARIATES}

if [ $? -ne 0 ]; then
    echo "ERROR: AEI-QTL Matrix Preparation failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Starting FastQTL for AEI-QTL Mapping ---"

# 2. Execution: FastQTL
${FASTQTL_BIN} \
    --vcf ${GENOTYPE_FILE} \
    --bed ${AEIQTL_PHENOTYPE} \
    --cov ${AEIQTL_COVARIATES} \
    --output ${FINAL_OUTPUT_FILE} \
    --window ${CIS_WINDOW_SIZE} \
    --permute ${PERMUTATIONS} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --chunk 1 1   # Run all analysis in one chunk

if [ $? -eq 0 ]; then
    echo "SUCCESS: Phase 7b AEI-QTL Mapping completed." | tee -a ${LOG_FILE}
    touch "${OUTPUT_DIR}/phase7b_success.flag"
else
    echo "ERROR: Phase 7b AEI-QTL Mapping failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Phase 7b Finished ---"#!/bin/bash
#
# SLURM wrapper for Phase 7b: AEI-QTL Association Mapping (FastQTL)
# Maps genetic variants to the Alu Editing Index (AEI).

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P7b_AEIQTL
#SBATCH --output=logs/phase7b_aeiqtl_%j.out
#SBATCH --error=logs/phase7b_aeiqtl_%j.err
#SBATCH --partition=bigmem        
#SBATCH --mem=64G                
#SBATCH --time=48:00:00          
#SBATCH --cpus-per-task=16        
#SBATCH --nodes=1

# --- Configuration: Paths ---
PREP_SCRIPT="./Python_scripts/prepare_aeiqtl_matrices.py"
FASTQTL_BIN="/path/to/FastQTL/executable"

# Input files (Outputs from Phase 6 and AEI_calculation)
INPUT_AEI_COVS="./phase6_normalized_edQTL/aei_covariate_matrix_p6a.tsv"
INPUT_COVS_DIR="./input_covariates"
INPUT_GENOTYPE_PCS="${INPUT_COVS_DIR}/genotype_pcs.tsv"
INPUT_PEER_FACTORS="${INPUT_COVS_DIR}/peer_factors_k60.tsv"

# Genotype and Output Paths
GENOTYPE_FILE="./genotype_data/genotype_all_individuals.vcf.gz"
OUTPUT_DIR="./phase7b_aeiqtl_results"
mkdir -p ${OUTPUT_DIR}

# AEI-QTL Specific Files
AEIQTL_PHENOTYPE="${OUTPUT_DIR}/aeiqtl_phenotype_p7b.tsv"
AEIQTL_COVARIATES="${OUTPUT_DIR}/aeiqtl_covariate_matrix_p7b.tsv"
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/aei_qtl_mapping_results_p7b.tsv.gz"
LOG_FILE="${OUTPUT_DIR}/aeiqtl_run_log.txt"

# Mapping Parameters (Same as P7)
CIS_WINDOW_SIZE=1000000
PERMUTATIONS=1000
# --- Setup ---
echo "--- Starting Phase 7b: AEI-QTL Mapping Preparation ---"
date | tee ${LOG_FILE}

# 1. Preparation (Run Python script)
python3 ${PREP_SCRIPT} \
    --input_aei_covariates ${INPUT_AEI_COVS} \
    --input_genotype_pcs ${INPUT_GENOTYPE_PCS} \
    --input_peer_factors ${INPUT_PEER_FACTORS} \
    --output_phenotype ${AEIQTL_PHENOTYPE} \
    --output_covariates ${AEIQTL_COVARIATES}

if [ $? -ne 0 ]; then
    echo "ERROR: AEI-QTL Matrix Preparation failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Starting FastQTL for AEI-QTL Mapping ---"

# 2. Execution: FastQTL
${FASTQTL_BIN} \
    --vcf ${GENOTYPE_FILE} \
    --bed ${AEIQTL_PHENOTYPE} \
    --cov ${AEIQTL_COVARIATES} \
    --output ${FINAL_OUTPUT_FILE} \
    --window ${CIS_WINDOW_SIZE} \
    --permute ${PERMUTATIONS} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --chunk 1 1   # Run all analysis in one chunk

if [ $? -eq 0 ]; then
    echo "SUCCESS: Phase 7b AEI-QTL Mapping completed." | tee -a ${LOG_FILE}
    touch "${OUTPUT_DIR}/phase7b_success.flag"
else
    echo "ERROR: Phase 7b AEI-QTL Mapping failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Phase 7b Finished ---"#!/bin/bash
#
# SLURM wrapper for Phase 7b: AEI-QTL Association Mapping (FastQTL)
# Maps genetic variants to the Alu Editing Index (AEI).

# --- SLURM Directives ---
#SBATCH --job-name=edQTL_P7b_AEIQTL
#SBATCH --output=logs/phase7b_aeiqtl_%j.out
#SBATCH --error=logs/phase7b_aeiqtl_%j.err
#SBATCH --partition=bigmem        
#SBATCH --mem=64G                
#SBATCH --time=48:00:00          
#SBATCH --cpus-per-task=16        
#SBATCH --nodes=1

# --- Configuration: Paths ---
PREP_SCRIPT="./Python_scripts/prepare_aeiqtl_matrices.py"
FASTQTL_BIN="/path/to/FastQTL/executable"

# Input files (Outputs from Phase 6 and AEI_calculation)
INPUT_AEI_COVS="./phase6_normalized_edQTL/aei_covariate_matrix_p6a.tsv"
INPUT_COVS_DIR="./input_covariates"
INPUT_GENOTYPE_PCS="${INPUT_COVS_DIR}/genotype_pcs.tsv"
INPUT_PEER_FACTORS="${INPUT_COVS_DIR}/peer_factors_k60.tsv"

# Genotype and Output Paths
GENOTYPE_FILE="./genotype_data/genotype_all_individuals.vcf.gz"
OUTPUT_DIR="./phase7b_aeiqtl_results"
mkdir -p ${OUTPUT_DIR}

# AEI-QTL Specific Files
AEIQTL_PHENOTYPE="${OUTPUT_DIR}/aeiqtl_phenotype_p7b.tsv"
AEIQTL_COVARIATES="${OUTPUT_DIR}/aeiqtl_covariate_matrix_p7b.tsv"
FINAL_OUTPUT_FILE="${OUTPUT_DIR}/aei_qtl_mapping_results_p7b.tsv.gz"
LOG_FILE="${OUTPUT_DIR}/aeiqtl_run_log.txt"

# Mapping Parameters (Same as P7)
CIS_WINDOW_SIZE=1000000
PERMUTATIONS=1000
# --- Setup ---
echo "--- Starting Phase 7b: AEI-QTL Mapping Preparation ---"
date | tee ${LOG_FILE}

# 1. Preparation (Run Python script)
python3 ${PREP_SCRIPT} \
    --input_aei_covariates ${INPUT_AEI_COVS} \
    --input_genotype_pcs ${INPUT_GENOTYPE_PCS} \
    --input_peer_factors ${INPUT_PEER_FACTORS} \
    --output_phenotype ${AEIQTL_PHENOTYPE} \
    --output_covariates ${AEIQTL_COVARIATES}

if [ $? -ne 0 ]; then
    echo "ERROR: AEI-QTL Matrix Preparation failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Starting FastQTL for AEI-QTL Mapping ---"

# 2. Execution: FastQTL
${FASTQTL_BIN} \
    --vcf ${GENOTYPE_FILE} \
    --bed ${AEIQTL_PHENOTYPE} \
    --cov ${AEIQTL_COVARIATES} \
    --output ${FINAL_OUTPUT_FILE} \
    --window ${CIS_WINDOW_SIZE} \
    --permute ${PERMUTATIONS} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --chunk 1 1   # Run all analysis in one chunk

if [ $? -eq 0 ]; then
    echo "SUCCESS: Phase 7b AEI-QTL Mapping completed." | tee -a ${LOG_FILE}
    touch "${OUTPUT_DIR}/phase7b_success.flag"
else
    echo "ERROR: Phase 7b AEI-QTL Mapping failed." | tee -a ${LOG_FILE}
    exit 1
fi

echo "--- Phase 7b Finished ---"
