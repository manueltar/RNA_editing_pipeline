#!/bin/bash
# Master SLURM Wrapper for RNA Editing Pipeline
# V17: Includes Robust Iteration, Error Checking, and Individual-Specific Logging/Cleanup

# --- Global Configuration ---
# 1. Define the path to your file containing one individual ID per line.
INDIVIDUALS_LIST="./individuals.txt"

# 2. Define the base path where the individual's cell-type BAMs are stored.
# ASSUMPTION: BAMs are located in $INPUT_DATA_ROOT_DIR/${INDIVIDUAL_ID}/
INPUT_DATA_ROOT_DIR="/path/to/your/input/bams/base_directory"

# 3. Define the base path for all pipeline output (intermediate and final)
OUTPUT_DATA_ROOT_DIR="/path/to/your/pipeline/output/base_directory"

# 4. Define the location of all phase scripts (Wrapper assumes they are here)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 5. Array Size Configuration (31 BAMs per individual)
ARRAY_SIZE=31

# 6. Pipeline Scripts - Ensure these names match the local files
SETUP_MASK_SCRIPT="${SCRIPT_DIR}/01_setup_references_mask_only_v2.sh"
PHASE2_CHUNK_PREP_SCRIPT="${SCRIPT_DIR}/02_index_pseudobamfiles_v3.sh"
PHASE2_REDML_SCRIPT="${SCRIPT_DIR}/call_red_ml_v2.sh"
PHASE2_REDITOOLS_SCRIPT="${SCRIPT_DIR}/call_reditools_v3.sh"
PHASE3_SCRIPT="${SCRIPT_DIR}/03_phase3_master_site_discovery_v2.sh"
PHASE4_SCRIPT="${SCRIPT_DIR}/04_per_cell_type_quantification_and_filter.sh"
PHASE5_SCRIPT="${SCRIPT_DIR}/05_final_aggregation.sh"
PHASE6_SCRIPT="${SCRIPT_DIR}/06_trait_normalization.sh"
CLEANUP_SCRIPT="${SCRIPT_DIR}/07_cleanup_intermediates.sh" # NEW CLEANUP SCRIPT

# --- Initial Setup ---
mkdir -p logs
echo "--- Starting RNA Editing Pipeline V17 (Processing Individuals from ${INDIVIDUALS_LIST}) ---"

# --- Main Iteration Loop ---
# This loop will run the entire pipeline sequence for each individual, one after the other.
while read INDIVIDUAL_ID; do
    if [[ -z "${INDIVIDUAL_ID}" ]]; then
        continue # Skip empty lines
    fi

    echo ""
    echo "=========================================================================="
    echo "STARTING PROCESSING FOR INDIVIDUAL: ${INDIVIDUAL_ID}"
    echo "=========================================================================="

    # --- Individual-Specific Path Setup ---
    INDIVIDUAL_LOG_DIR="logs/${INDIVIDUAL_ID}"
    INDIVIDUAL_OUTPUT_DIR="${OUTPUT_DATA_ROOT_DIR}/${INDIVIDUAL_ID}"
    INPUT_BAM_ROOT="${INPUT_DATA_ROOT_DIR}/${INDIVIDUAL_ID}/cell_type_bams"
    
    mkdir -p "${INDIVIDUAL_LOG_DIR}"
    mkdir -p "${INDIVIDUAL_OUTPUT_DIR}"

    # --- Function to check sbatch submission success ---
    sbatch_check () {
        local job_id="$1"
        local job_name="$2"
        if [[ -z "${job_id}" ]]; then
            echo "ERROR: SLURM failed to submit ${job_name}. Pipeline terminated for ${INDIVIDUAL_ID}."
            exit 1
        fi
        echo "Submitted ${job_name} Job ID: ${job_id}"
    }

    # --------------------------------------------------------------------------
    # PHASE 1: SETUP (SERIAL)
    # --------------------------------------------------------------------------
    echo "Submitting Phase 1: Reference Mask Setup (Serial Job)"
    JOB_ID_P1_SETUP=$(sbatch \
        --job-name=P1_Setup_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p1_%j.out \
        ${SETUP_MASK_SCRIPT} ${INDIVIDUAL_ID} ${INPUT_BAM_ROOT} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P1_SETUP}" "Phase 1 Setup"


    # --------------------------------------------------------------------------
    # PHASE 2: DATA PREP & CALLING (PARALLEL ARRAY JOBS)
    # --------------------------------------------------------------------------

    # P2: Chunk Preparation (N=31 separate BAM files)
    echo "Submitting Phase 2: Chunk Preparation (${ARRAY_SIZE} Array Tasks)"
    JOB_ID_P2_CHUNK_PREP=$(sbatch \
        --dependency=afterok:${JOB_ID_P1_SETUP} \
        --job-name=P2_ChunkPrep_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p2_chunk_%A_%a.out \
        --array=1-${ARRAY_SIZE}%10 \
        ${PHASE2_CHUNK_PREP_SCRIPT} ${INDIVIDUAL_ID} ${INPUT_BAM_ROOT} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P2_CHUNK_PREP}" "Phase 2 Chunk Prep Array"

    # P2: RNA Editing Calling (REDML)
    echo "Submitting Phase 2: RED-ML Calling (${ARRAY_SIZE} Array Tasks)"
    JOB_ID_P2_REDML=$(sbatch \
        --dependency=afterok:${JOB_ID_P2_CHUNK_PREP}_${ARRAY_SIZE} \
        --job-name=P2_REDML_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p2_redml_%A_%a.out \
        --array=1-${ARRAY_SIZE} \
        ${PHASE2_REDML_SCRIPT} ${INDIVIDUAL_ID} ${INPUT_BAM_ROOT} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P2_REDML}" "Phase 2 RED-ML Array"

    # P2: RNA Editing Calling (REDItools)
    echo "Submitting Phase 2: REDItools Calling (${ARRAY_SIZE} Array Tasks)"
    JOB_ID_P2_REDITOOLS=$(sbatch \
        --dependency=afterok:${JOB_ID_P2_CHUNK_PREP}_${ARRAY_SIZE} \
        --job-name=P2_REDItools_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p2_reditools_%A_%a.out \
        --array=1-${ARRAY_SIZE} \
        ${PHASE2_REDITOOLS_SCRIPT} ${INDIVIDUAL_ID} ${INPUT_BAM_ROOT} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P2_REDITOOLS}" "Phase 2 REDItools Array"

    # --------------------------------------------------------------------------
    # PHASE 3: MASTER SITE DISCOVERY (SERIAL)
    # --------------------------------------------------------------------------
    echo "Submitting Phase 3: Master Site Discovery (Serial Job)"
    JOB_ID_P3=$(sbatch \
        --dependency=afterok:${JOB_ID_P2_REDML}_${ARRAY_SIZE}:${JOB_ID_P2_REDITOOLS}_${ARRAY_SIZE} \
        --job-name=P3_Master_Sites_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p3_%j.out \
        ${PHASE3_SCRIPT} ${INDIVIDUAL_ID} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P3}" "Phase 3 Master Sites"

    # --------------------------------------------------------------------------
    # PHASE 4: PER-CELL QUANTIFICATION & FILTERING (SERIAL)
    # --------------------------------------------------------------------------
    echo "Submitting Phase 4: Per-Cell Quantification & Final Filtering (Serial Job)"
    JOB_ID_P4=$(sbatch \
        --dependency=afterok:${JOB_ID_P3} \
        --job-name=P4_Quant_Filter_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p4_%j.out \
        ${PHASE4_SCRIPT} ${INDIVIDUAL_ID} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P4}" "Phase 4 Quantification"

    # --------------------------------------------------------------------------
    # PHASE 5: FINAL AGGREGATION (SERIAL)
    # --------------------------------------------------------------------------
    echo "Submitting Phase 5: Final Aggregation (Serial Job)"
    JOB_ID_P5=$(sbatch \
        --dependency=afterok:${JOB_ID_P4} \
        --job-name=P5_Aggregation_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p5_%j.out \
        ${PHASE5_SCRIPT} ${INDIVIDUAL_ID} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P5}" "Phase 5 Aggregation"

    # --------------------------------------------------------------------------
    # PHASE 6: TRAIT NORMALIZATION (SERIAL - NEW)
    # --------------------------------------------------------------------------
    echo "Submitting Phase 6: Trait Normalization (Serial Job)"
    JOB_ID_P6=$(sbatch \
        --dependency=afterok:${JOB_ID_P5} \
        --job-name=P6_Normalization_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p6_%j.out \
        ${PHASE6_SCRIPT} ${INDIVIDUAL_ID} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P6}" "Phase 6 Normalization"

    # --------------------------------------------------------------------------
    # PHASE 7: CLEANUP (SERIAL - NEW)
    # --------------------------------------------------------------------------
    # This step ensures intermediate files are removed for the current individual
    # before the next individual begins, managing HPC storage.
    echo "Submitting Phase 7: Cleanup (Serial Job)"
    JOB_ID_P7_CLEANUP=$(sbatch \
        --dependency=afterok:${JOB_ID_P6} \
        --job-name=P7_Cleanup_${INDIVIDUAL_ID} \
        --output=${INDIVIDUAL_LOG_DIR}/p7_%j.out \
        ${CLEANUP_SCRIPT} ${INDIVIDUAL_ID} ${INDIVIDUAL_OUTPUT_DIR} | awk '{print $4}')
    sbatch_check "${JOB_ID_P7_CLEANUP}" "Phase 7 Cleanup"

    echo "--- Individual ${INDIVIDUAL_ID} job chain submitted. Final job ID: ${JOB_ID_P7_CLEANUP} ---"

    # CRITICAL: Wait for the P7 cleanup job to complete before starting the next individual's chain.
    # This enforces the required sequential processing of individuals.
    srun --dependency=afterok:${JOB_ID_P7_CLEANUP} --job-name=Barrier_For_${INDIVIDUAL_ID} --ntasks=1 --time=00:01:00 echo "Barrier cleared for next individual."
    
done < "${INDIVIDUALS_LIST}"

echo "=========================================================================="
echo "All individuals submitted for processing. Check logs directory for output."
echo "=========================================================================="
