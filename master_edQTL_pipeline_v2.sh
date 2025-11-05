#!/bin/bash
#
# MASTER SLURM WRAPPER: edQTL Statistical Pipeline (Phases 5 through 8)
#
# This script submits the sequential jobs necessary for feature selection, 
# AEI covariate calculation, normalization, association mapping, and multiple testing correction.

# --- Configuration: Script Paths ---
PHASE5_SCRIPT="./master_sbatch_scripts/run_phase5_collation.sh"
AEI_SCRIPT="./master_sbatch_scripts/run_AEI_calculation_array.sh" 
PHASE6_SCRIPT="./master_sbatch_scripts/run_phase6_processing.sh"
PHASE7_SCRIPT="./master_sbatch_scripts/run_phase7_edqtl_mapping.sh"
PHASE8_SCRIPT="./master_sbatch_scripts/run_phase8_qvalue_filter.sh"

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
# This phase calculates the AEI matrix which is required by Phase 6.
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
# --- 4. Phase 7: edQTL Association Mapping (FastQTL) ---
# Dependency: Requires Phase 6 (INT-transformed phenotype and merged covariates).
# ----------------------------------------------------------------------
echo "Submitting Phase 7: FastQTL Mapping, dependent on Phase 6 success..."
PHASE7_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE6_JOB_ID} ${PHASE7_SCRIPT})

if [ -z "${PHASE7_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 7. Aborting."
    exit 1
fi
echo "Phase 7 submitted (Job ID: ${PHASE7_JOB_ID}, Dependent on ${PHASE6_JOB_ID})"

# ----------------------------------------------------------------------
# --- 5. Phase 8: Multiple Testing Correction (Q-value Filtering) ---
# Dependency: Requires Phase 7 (raw FastQTL permutation results).
# ----------------------------------------------------------------------
echo "Submitting Phase 8: FDR Correction and Lead SNP Identification, dependent on Phase 7 success..."
PHASE8_JOB_ID=$(sbatch --parsable --dependency=afterok:${PHASE7_JOB_ID} ${PHASE8_SCRIPT})

if [ -z "${PHASE8_JOB_ID}" ]; then
    echo "ERROR: Failed to submit Phase 8. Aborting."
    exit 1
fi
echo "Phase 8 submitted (Job ID: ${PHASE8_JOB_ID}, Dependent on ${PHASE7_JOB_ID})"

echo "--- All statistical pipeline jobs submitted successfully ---"
echo "Final job chain: P5 (${PHASE5_JOB_ID}) & AEI (${AEI_JOB_ID}) -> P6 (${PHASE6_JOB_ID}) -> P7 (${PHASE7_JOB_ID}) -> P8 (${PHASE8_JOB_ID})"
