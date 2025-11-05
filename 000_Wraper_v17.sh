#!/bin/bash
# Master SLURM Wrapper for RNA Editing Pipeline 
# V5: Includes Reworked Downstream Steps (P3/P4), P6 Trait Normalization (via downstream master script),
#     AND a final **Cleanup Phase** to manage file quota after P4.

# --- Define Paths (Assuming your scripts are in the same directory) ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Phase 1 Scripts
SETUP_MASK_SCRIPT="${SCRIPT_DIR}/01_setup_references_mask_only_v2.sh"

# Phase 2 Scripts (Chunk Prep)
PHASE2_CHUNK_PREP_SCRIPT="${SCRIPT_DIR}/02_index_pseudobamfiles_v2.sh"

# Phase 2 Scripts (RNA Editing Calling)
PHASE2_REDML_SCRIPT="${SCRIPT_DIR}/call_red_ml_v3.sh"
PHASE2_REDITOOLS_SCRIPT="${SCRIPT_DIR}/call_reditools_v5.sh"

# Phase 3 Scripts (Master Site Discovery - REWORKED)
PHASE3_SCRIPT="${SCRIPT_DIR}/03_phase3_master_site_discovery_v5.sh"

# Phase 4 Scripts (Per-Cell Quantification & Filter - REWORKED)
PHASE4_SCRIPT="${SCRIPT_DIR}/04_per_cell_type_quantification_and_filter_v3.sh"

# --- NEW Cleanup Script ---
CLEANUP_SCRIPT="${SCRIPT_DIR}/run_cleanup_p1_p4.sh" 

# --- Array Size Configuration (31 BAMs per individual) ---
ARRAY_SIZE=31

# --- Log Directory Setup ---
mkdir -p logs

echo "--- Starting RNA Editing Pipeline V5 (Includes P1-P4 Cleanup) ---"
echo "Processing with Array Size: ${ARRAY_SIZE}"

# --------------------------------------------------------------------------
# PHASE 1: SETUP
# --------------------------------------------------------------------------

# P1: Setup References and Blacklist (Serial Job)
echo "Submitting Phase 1: Reference Mask Setup (Serial Job)"
JOB_ID_P1_SETUP=$(sbatch ${SETUP_MASK_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 1 Setup Job ID: ${JOB_ID_P1_SETUP}"

# --------------------------------------------------------------------------
# PHASE 2: DATA PREP & CALLING (Parallel Array Jobs)
# --------------------------------------------------------------------------

# P2: Chunk Preparation (N=${ARRAY_SIZE} separate BAM files)
echo "Submitting Phase 2: Chunk Preparation (${ARRAY_SIZE} Array Tasks)"
JOB_ID_P2_CHUNK_PREP=$(sbatch \
    --dependency=afterok:${JOB_ID_P1_SETUP} \
    --job-name=P2_ChunkPrep \
    --array=1-${ARRAY_SIZE}%10 \
    ${PHASE2_CHUNK_PREP_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 2 Chunk Preparation Array Job ID: ${JOB_ID_P2_CHUNK_PREP}"

# P2: RNA Editing Calling (REDML)
echo "Submitting Phase 2: RED-ML Calling (${ARRAY_SIZE} Array Tasks)"
JOB_ID_P2_REDML=$(sbatch \
    --dependency=afterok:${JOB_ID_P2_CHUNK_PREP} \
    --job-name=P2_REDML_Call \
    --array=1-${ARRAY_SIZE} \
    ${PHASE2_REDML_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 2 RED-ML Array Job ID: ${JOB_ID_P2_REDML}"

# P2: RNA Editing Calling (REDItools)
echo "Submitting Phase 2: REDItools Calling (${ARRAY_SIZE} Array Tasks)"
JOB_ID_P2_REDITOOLS=$(sbatch \
    --dependency=afterok:${JOB_ID_P2_CHUNK_PREP} \
    --job-name=P2_REDItools_Call \
    --array=1-${ARRAY_SIZE} \
    ${PHASE2_REDITOOLS_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 2 REDItools Array Job ID: ${JOB_ID_P2_REDITOOLS}"

# --------------------------------------------------------------------------
# PHASE 3: MASTER SITE DISCOVERY (Serial Job)
# --------------------------------------------------------------------------

# P3: Aggregates all raw chunk files to establish a universal list of consensus sites.
# Dependency: Wait for ALL tasks from P2_REDML and P2_REDItools
echo "Submitting Phase 3: Master Site Discovery (Serial Job)"
JOB_ID_P3=$(sbatch \
    --dependency=afterok:${JOB_ID_P2_REDML}:${JOB_ID_P2_REDITOOLS} \
    --job-name=P3_Master_Sites \
    ${PHASE3_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 3 Job ID: ${JOB_ID_P3}"

# --------------------------------------------------------------------------
# PHASE 4: PER-CELL QUANTIFICATION, FILTERING, AND ANNOTATION (Serial Job)
# --------------------------------------------------------------------------

# P4: Takes the master site list and queries each of the 31 BAMs for final metrics and filtering.
# Dependency: Phase 3
echo "Submitting Phase 4: Per-Cell Quantification & Final Filtering (Serial Job)"
JOB_ID_P4=$(sbatch \
    --dependency=afterok:${JOB_ID_P3} \
    --job-name=P4_Quant_Filter \
    ${PHASE4_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 4 Job ID: ${JOB_ID_P4}"

# --------------------------------------------------------------------------
# NEW CLEANUP PHASE (Serial Job)
# --------------------------------------------------------------------------

# Cleanup: Erases intermediate files to free up disk space.
# Dependency: Requires Phase 4 to complete successfully.
echo "Submitting Cleanup Job, dependent on Phase 4 success..."
JOB_ID_CLEANUP=$(sbatch \
    --dependency=afterok:${JOB_ID_P4} \
    --job-name=P1_P4_Cleanup \
    ${CLEANUP_SCRIPT} | awk '{print $4}')
echo "Submitted Cleanup Job ID: ${JOB_ID_CLEANUP}"


echo "--- All P1-P4 Processing Jobs Submitted. ---"
echo "Job Chain: P1 -> P2_ChunkPrep -> (P2_REDML & P2_REDItools) -> P3 -> P4 -> Cleanup (${JOB_ID_CLEANUP})"
