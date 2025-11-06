#!/bin/bash
#
# MASTER SLURM WRAPPER: edQTL Statistical Pipeline (Phases 5 through 8)
# V3: Adds Phase 7b for AEI-QTL mapping in parallel with Phase 7.
#     Phase 8 is updated to depend on both P7 and P7b for combined FDR correction.

# --- Configuration: Script Paths ---
PHASE5_SCRIPT="./master_sbatch_scripts/run_phase5_collation_v2.sh"
AEI_SCRIPT="./master_sbatch_scripts/run_AEI_calculation_array.sh" 
PHASE6_SCRIPT="./master_sbatch_scripts/run_phase6_processing_v2.sh"
PHASE7_SCRIPT="./master_sbatch_scripts/run_phase7_edqtl_mapping_v2.sh"
PHASE7B_SCRIPT="./master_sbatch_scripts/run_phase7b_aeiqtl_mapping_v3.sh" 
PHASE8_SCRIPT="./master_sbatch_scripts/run_phase8_qvalue_filter_v2.sh"

# --- Setup ---
echo "--- Starting Master edQTL Statistical Pipeline Submission ---"

# ----------------------------------------------------------------------
# --- 1. Phase 5: Feature Selection and Collation ---
# ----------------------------------------------------------------------
echo "Submitting Phase 5: Collation and Feature Selection..."
PHASE5_JOB_ID=$(sbatch --parsable ${PHASE5_SCRIPT})

if [ -z "${PHASE5_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 5. Aborting."
    exit 1
fi
echo "Phase 5 submitted (Job ID: ${PHASE5_JOB_ID})"

# ----------------------------------------------------------------------
# --- 2. AEI_calculation: Alu Editing Index (Covariate Generation) ---
# Dependency: Runs after Phase 5 completes (proxy for pipeline stability).
# ----------------------------------------------------------------------
echo "Submitting AEI Calculation and Collation, dependent on Phase 5 success..."
AEI_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE5_JOB_ID} ${AEI_SCRIPT})

if [ -z "${AEI_JOB_ID}" ]; then
    echo "ERROR: Failed to submit AEI Calculation. Aborting."
    exit 1
fi
echo "AEI Calculation submitted (Job ID: ${AEI_JOB_ID}, Dependent on ${PHASE5_JOB_ID})"

# ----------------------------------------------------------------------
# --- 3. Phase 6: Normalization and Covariate Merge ---
# Dependency: Requires BOTH Phase 5 (feature data) and AEI_calculation (AEI covariate data).
# ----------------------------------------------------------------------
echo "Submitting Phase 6: Normalization and Covariate Merge, dependent on Phase 5 AND AEI success..."
PHASE6_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE5_JOB_ID}:${AEI_JOB_ID} ${PHASE6_SCRIPT})

if [ -z "${PHASE6_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 6. Aborting."
    exit 1
fi
echo "Phase 6 submitted (Job ID: ${PHASE6_JOB_ID}, Dependent on ${PHASE5_JOB_ID} and ${AEI_JOB_ID})"

# ----------------------------------------------------------------------
# --- 4. Phase 7 (edQTL) and Phase 7b (AEI-QTL) - PARALLEL MAPPING ---
# Dependency: Both require Phase 6 to produce the input matrices.
# ----------------------------------------------------------------------
echo "Submitting Phase 7: edQTL Mapping, dependent on Phase 6 success..."
PHASE7_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE6_JOB_ID} ${PHASE7_SCRIPT})

if [ -z "${PHASE7_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 7. Aborting."
    exit 1
fi
echo "Phase 7 submitted (Job ID: ${PHASE7_JOB_ID}, Dependent on ${PHASE6_JOB_ID})"

echo "Submitting Phase 7b: AEI-QTL Mapping, dependent on Phase 6 success..."
PHASE7B_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE6_JOB_ID} ${PHASE7B_SCRIPT})

if [ -z "${PHASE7B_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 7b. Aborting."
    exit 1
fi
echo "Phase 7b submitted (Job ID: ${PHASE7B_JOB_ID}, Dependent on ${PHASE6_JOB_ID})"

# ----------------------------------------------------------------------
# --- 5. Phase 8: Multiple Testing Correction (Q-value Filtering) ---
# Dependency: Requires BOTH Phase 7 (edQTL) and Phase 7b (AEI-QTL) to finish.
# ----------------------------------------------------------------------
echo "Submitting Phase 8: FDR Correction and Lead SNP Identification, dependent on Phase 7 AND 7b success..."
PHASE8_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE7_JOB_ID}:${PHASE7B_JOB_ID} ${PHASE8_SCRIPT})

if [ -z "${PHASE8_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 8. Aborting."
    exit 1
fi
echo "Phase 8 submitted (Job ID: ${PHASE8_JOB_ID}, Dependent on ${PHASE7_JOB_ID} and ${PHASE7B_JOB_ID})"

echo "--- All statistical pipeline jobs submitted successfully ---"
echo "Final job chain: P5 & AEI -> P6 -> P7 & P7b (Parallel) -> P8"
