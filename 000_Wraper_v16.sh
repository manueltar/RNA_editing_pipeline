#!/bin/bash
# Master SLURM Wrapper for RNA Editing Pipeline 
# V5: Includes Reworked Downstream Steps (P3/P4) and NEW P6 Trait Normalization.

# --- Define Paths (Assuming your scripts are in the same directory) ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Phase 1 Scripts
SETUP_MASK_SCRIPT="${SCRIPT_DIR}/01_setup_references_mask_only_v2.sh"

# Phase 2 Scripts (Chunk Prep)
PHASE2_CHUNK_PREP_SCRIPT="${SCRIPT_DIR}/02_index_pseudobamfiles_v2.sh"

# Phase 2 Scripts (RNA Editing Calling)
PHASE2_REDML_SCRIPT="${SCRIPT_DIR}/call_red_ml_v2.sh"
PHASE2_REDITOOLS_SCRIPT="${SCRIPT_DIR}/call_reditools_v3.sh"

# Phase 3 Scripts (Master Site Discovery - REWORKED)
PHASE3_SCRIPT="${SCRIPT_DIR}/03_phase3_master_site_discovery_v2.sh"

# Phase 4 Scripts (Per-Cell Quantification & Filter - REWORKED)
PHASE4_SCRIPT="${SCRIPT_DIR}/04_per_cell_type_quantification_and_filter.sh"

# Phase 5 Scripts (Final Aggregation)
PHASE5_SCRIPT="${SCRIPT_DIR}/05_final_aggregation.sh" # <--- UPDATED SCRIPT NAME

# Phase 6 Scripts (Trait Normalization - NEW)
PHASE6_SCRIPT="${SCRIPT_DIR}/06_trait_normalization.sh" # <--- NEW SCRIPT

# --- Array Size Configuration (31 BAMs per individual) ---
ARRAY_SIZE=31

# --- Log Directory Setup ---
mkdir -p logs

echo "--- Starting RNA Editing Pipeline V5 (Includes P6 Trait Normalization) ---"
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
    --dependency=afterok:${JOB_ID_P2_CHUNK_PREP}_${ARRAY_SIZE} \
    --job-name=P2_REDML_Call \
    --array=1-${ARRAY_SIZE} \
    ${PHASE2_REDML_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 2 RED-ML Array Job ID: ${JOB_ID_P2_REDML}"

# P2: RNA Editing Calling (REDItools)
echo "Submitting Phase 2: REDItools Calling (${ARRAY_SIZE} Array Tasks)"
JOB_ID_P2_REDITOOLS=$(sbatch \
    --dependency=afterok:${JOB_ID_P2_CHUNK_PREP}_${ARRAY_SIZE} \
    --job-name=P2_REDItools_Call \
    --array=1-${ARRAY_SIZE} \
    ${PHASE2_REDITOOLS_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 2 REDItools Array Job ID: ${JOB_ID_P2_REDITOOLS}"

# --------------------------------------------------------------------------
# PHASE 3: MASTER SITE DISCOVERY (Serial Job)
# --------------------------------------------------------------------------

# P3: Aggregates all raw chunk files to establish a universal list of consensus sites.
echo "Submitting Phase 3: Master Site Discovery (Serial Job)"
JOB_ID_P3=$(sbatch \
    --dependency=afterok:${JOB_ID_P2_REDML}_${ARRAY_SIZE}:${JOB_ID_P2_REDITOOLS}_${ARRAY_SIZE} \
    --job-name=P3_Master_Sites \
    ${PHASE3_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 3 Job ID: ${JOB_ID_P3}"

# --------------------------------------------------------------------------
# PHASE 4: PER-CELL QUANTIFICATION, FILTERING, AND ANNOTATION (Serial Job)
# --------------------------------------------------------------------------

# P4: Takes the master site list and queries each of the 31 BAMs for final metrics and filtering.
echo "Submitting Phase 4: Per-Cell Quantification & Final Filtering (Serial Job)"
JOB_ID_P4=$(sbatch \
    --dependency=afterok:${JOB_ID_P3} \
    --job-name=P4_Quant_Filter \
    ${PHASE4_SCRIPT} | awk '{print $4}')
echo "Submitted Phase 4 Job ID: ${JOB_ID_P4}"

# --------------------------------------------------------------------------
# PHASE 5: FINAL AGGREGATION (Serial Job)
# --------------------------------------------------------------------------

# P5: Aggregates the P4 quantification matrix into gene-level AvgER traits.
echo "Submitting Phase 5: Final Aggregation (Serial Job)"
JOB_ID_P5=$(sbatch \
    --dependency=afterok:${JOB_ID_P4} \
    --job-name=P5_Aggregation \
    ${PHASE5_SCRIPT} | awk '{print $4}') # <--- SCRIPT NAME CORRECTED HERE
echo "Submitted Phase 5 Job ID: ${JOB_ID_P5}"

# --------------------------------------------------------------------------
# PHASE 6: TRAIT NORMALIZATION (Serial Job - NEW)
# --------------------------------------------------------------------------

# P6: Applies Inverse Normal Transformation (INT) to the gene-level AvgER traits.
echo "Submitting Phase 6: Trait Normalization (Serial Job)"
JOB_ID_P6=$(sbatch \
    --dependency=afterok:${JOB_ID_P5} \
    --job-name=P6_Normalization \
    ${PHASE6_SCRIPT} | awk '{print $4}') # <--- NEW DEPENDENT JOB
echo "Submitted Phase 6 Job ID: ${JOB_ID_P6}"

echo "--- All jobs submitted. Final job ID: ${JOB_ID_P6} ---"
