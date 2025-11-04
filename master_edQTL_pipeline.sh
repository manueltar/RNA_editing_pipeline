#!/bin/bash
# --- Top-level Master Script for edQTL Pipeline ---
# This script manages the flow between major phases.

# Directories (Consistent with Phase 5 script)
PHASE5_OUTPUT_DIR="./phase5_edQTL_features"
PHASE5_FLAG="${PHASE5_OUTPUT_DIR}/phase5_success.flag"
PHASE6_SCRIPT="./run_phase6_INT_covariates.sh" # Placeholder for the next phase

# --- Phase 5: Population Feature Selection (Collation) ---

echo "--- Submitting Phase 5: edQTL Feature Collation Job ---"

# Submit the heavy collation job
PHASE5_JOB_ID=$(sbatch run_phase5_collation.sh | awk '{print $4}')

if [ -z "$PHASE5_JOB_ID" ]; then
    echo "FATAL ERROR: Failed to submit Phase 5 job. Exiting."
    exit 1
fi

echo "Phase 5 Job ID: $PHASE5_JOB_ID"

# Wait for Phase 5 to complete successfully
echo "Waiting for Phase 5 job ($PHASE5_JOB_ID) to complete successfully..."
# This command waits for the dependency job to finish
squeue -j $PHASE5_JOB_ID -o %T | grep -q "RUNNING\|PENDING"
while [ $? -eq 0 ]; do
    sleep 300 # Wait 5 minutes between checks
    squeue -j $PHASE5_JOB_ID -o %T | grep -q "RUNNING\|PENDING"
done

# Check the success flag created by the Phase 5 script
if [ -f "$PHASE5_FLAG" ]; then
    echo "SUCCESS: Phase 5 edQTL Feature Selection completed. Proceeding to Phase 6 setup."
    
    # --- Phase 6: Next Step (e.g., INT, Covariate Prep, and FastQTL setup) ---
    
    # Note: Phase 6 would typically be submitted here, dependent on the completion of Phase 5.
    # We will generate this script later.
    # sbatch --depend=afterok:$PHASE5_JOB_ID $PHASE6_SCRIPT

else
    echo "ERROR: Phase 5 failed. The success flag file ($PHASE5_FLAG) was not found. Please check SLURM logs."
fi

echo "--- Master Script Finished ---"
