#!/bin/bash
# 
# Master Pipeline Script for edQTL Analysis (Phases 5 & 6)
# This script orchestrates the submission of sequential SLURM jobs, ensuring
# that the feature selection (P5) completes before normalization (P6) begins.

# --- Configuration ---
PHASE5_SCRIPT="./run_phase5_collation_v2.sh"
PHASE6_SCRIPT="./run_phase6_processing.sh"

PHASE5_FLAG="./phase5_edQTL_features/phase5_success.flag"
PHASE6_FLAG="./phase6_normalized_edQTL/phase6_success.flag"

# --- Setup ---
echo "--- Starting Master edQTL Pipeline ---"

# --- Phase 5: Collation and Feature Selection ---
if [ -f "$PHASE5_FLAG" ]; then
    echo "Phase 5 success flag found. Skipping submission."
else
    echo "Submitting Phase 5: edQTL Feature Collation and Selection..."
    # The submission uses sbatch and captures the Job ID
    JOBID_P5=$(sbatch --parsable ${PHASE5_SCRIPT})
    if [ -z "$JOBID_P5" ]; then
        echo "ERROR: Phase 5 submission failed."
        exit 1
    fi
    echo "Phase 5 submitted with Job ID: $JOBID_P5. Waiting for completion..."

    # Use a loop to check for the success flag
    while [ ! -f "$PHASE5_FLAG" ]; do
        # Check job status to ensure it hasn't failed or been canceled
        if ! squeue -j "$JOBID_P5" &> /dev/null; then
            echo "ERROR: Phase 5 job $JOBID_P5 failed or completed without success flag. Check logs/ for details."
            exit 1
        fi
        sleep 60 # Wait for 60 seconds before checking again
    done
    echo "Phase 5 completed successfully."
fi

# --- Phase 6: Normalization and Covariate Adjustment ---
if [ -f "$PHASE6_FLAG" ]; then
    echo "Phase 6 success flag found. Skipping submission."
else
    echo "Submitting Phase 6: Inverse Normal Transformation and Covariate Prep..."
    
    # Submit Phase 6, ensuring it runs only after Phase 5 has finished.
    # Note: If Phase 5 was skipped (flag found), Phase 6 runs immediately.
    JOBID_P6=$(sbatch --parsable --dependency=afterok:$JOBID_P5 ${PHASE6_SCRIPT})
    if [ -z "$JOBID_P6" ]; then
        echo "ERROR: Phase 6 submission failed."
        exit 1
    fi
    echo "Phase 6 submitted with Job ID: $JOBID_P6. Waiting for completion..."

    # Use a loop to check for the success flag
    while [ ! -f "$PHASE6_FLAG" ]; do
        # Check job status
        if ! squeue -j "$JOBID_P6" &> /dev/null; then
            echo "ERROR: Phase 6 job $JOBID_P6 failed or completed without success flag. Check logs/ for details."
            exit 1
        fi
        sleep 60 # Wait for 60 seconds before checking again
    done
    echo "Phase 6 completed successfully. Files are ready for FastQTL."
fi

echo "--- Master edQTL Pipeline Finished ---"
