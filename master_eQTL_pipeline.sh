#!/bin/bash
# Master script to control the execution flow of the eQTL pipeline.
# Assumes Phases 1-4 (per-individual processing) are complete for all 6000 samples.
 
PHASE5_SCRIPT="./run_phase5_collation.sh"
PHASE5_OUTPUT_FLAG="./phase5_eQTL_features/phase5_success.flag"
PHASE6_SCRIPT="./run_phase6_INT_and_covariates.sh" # Placeholder for the next step
PHASE7_SCRIPT="./run_phase7_fastQTL.sh"             # Placeholder for the final eQTL run

echo "--- Starting Master eQTL Pipeline Control ---"
date

# --- PHASE 5: Population Collation and Feature Selection ---
echo -e "\n[PHASE 5] Submitting Population Collation Job..."

# Check if Phase 5 has already completed successfully
if [ -f "$PHASE5_OUTPUT_FLAG" ]; then
    echo "Phase 5 Collation already marked as SUCCESS. Skipping submission."
else
    # Submit the Phase 5 job
    PHASE5_JOB_ID=$(sbatch --parsable $PHASE5_SCRIPT)
    if [ $? -eq 0 ]; then
        echo "Phase 5 job submitted successfully. JOB ID: $PHASE5_JOB_ID"
        
        # We need to wait for Phase 5 to finish before continuing
        echo "Waiting for Phase 5 job ($PHASE5_JOB_ID) to complete..."
        # You can use 'scontrol' or a loop to monitor status, or submit Phase 6 with a dependency
        
        # Option 1: Submit Phase 6 with a dependency (Recommended for large pipelines)
        # We will use this dependency to ensure serial execution
        echo "Phase 6 submission will be dependent on Phase 5 completion..."
        
        # --- PHASE 6: INT Transformation and Covariate Preparation (Placeholder) ---
        # Note: Replace --dependency with the actual script name/logic when ready for Phase 6
        
        # PHASE6_JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE5_JOB_ID $PHASE6_SCRIPT)
        # echo "Phase 6 job submitted successfully. JOB ID: $PHASE6_JOB_ID"
        
    else
        echo "ERROR: Phase 5 job submission failed."
        exit 1
    fi
fi

echo -e "\nMaster script finished submitting initial jobs."
echo "Please monitor job status on your cluster."
