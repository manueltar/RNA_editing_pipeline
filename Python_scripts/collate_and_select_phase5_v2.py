#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
import sys
import os

# This script implements the definitive Li et al. (2022) edQTL feature selection strategy 
# by collating data across ALL individuals and then selecting the most active site 
# per (Gene, CellType) feature based on the population median of raw editing ratios.

# --- Main Function ---
def run_phase5_collation_and_selection(args):
    """
    1. Loads all individual Phase 4 matrices.
    2. Collates them into a single population-level table (long format).
    3. **APPLIES SAMPLE SIZE FILTER (N >= args.min_samples)**.
    4. Performs the Representative Site selection (highest median raw ER) 
       for each unique (Gene, CellType) feature across the entire population.
    """
    MIN_SAMPLES = args.min_samples # Now read from argparse
    
    input_pattern = os.path.join(args.input_dir, args.file_pattern)
    input_files = glob.glob(input_pattern)
    
    if not input_files:
        print(f"FATAL ERROR: No files found matching pattern: {input_pattern}", file=sys.stderr)
        return

    print(f"Found {len(input_files)} individual files for collation.")
    all_data = []

    # 1. Collation (Iterate and Stack) - Converting from wide (by CellType) to long (by Sample/Individual)
    for i, file_path in enumerate(input_files):
        if (i + 1) % 500 == 0 or i == 0 or i == len(input_files) - 1:
            print(f"  Processing file {i+1}/{len(input_files)}: {os.path.basename(file_path)}")
        try:
            # Read the Phase 4 file
            df = pd.read_csv(file_path, sep='\t', index_col='SiteID', na_values=['NA'])
            
            # --- Essential: Filter for Sites that Passed Global QC ---
            df_filtered = df[df['GlobalFilterStatus'] == 'PASS'].copy()
            
            # Extract individual ID from the filename (e.g., IID_final_editing_matrix_p4.tsv -> IID)
            individual_id = os.path.basename(file_path).replace(args.file_pattern.replace('*', ''), '')
            
            # Identify relevant columns (Editing Ratios and Gene annotation)
            er_cols = [col for col in df_filtered.columns if col.endswith('_ER')]
            metadata_cols = ['Phase3_Gene']
            
            # Convert the individual's wide table to long format
            df_long = df_filtered[metadata_cols + er_cols].reset_index().melt(
                id_vars=['SiteID', 'Phase3_Gene'], 
                value_vars=er_cols, 
                var_name='CellType_ER', 
                value_name=individual_id
            )
            
            # Clean up and prepare the long format for merging
            df_long['CellType'] = df_long['CellType_ER'].str.replace('_ER', '')
            df_long = df_long.drop(columns=['CellType_ER'])
            
            # Set index for robust merging: (SiteID, Gene, CellType)
            df_long.set_index(['SiteID', 'Phase3_Gene', 'CellType'], inplace=True)
            all_data.append(df_long[[individual_id]])

        except Exception as e:
            print(f"WARNING: Skipping file {file_path} due to error: {e}", file=sys.stderr)

    if not all_data:
        print("No data successfully processed. Exiting.", file=sys.stderr)
        return

    # Combine all long dataframes into one massive population matrix
    population_matrix_raw = pd.concat(all_data, axis=1)
    
    print(f"\n--- Collation Complete ---")
    print(f"Raw population matrix shape: {population_matrix_raw.shape}")
    
    # 2. **Sample Size Filter (NEW STEP)**
    print(f"Applying Sample Size Filter: Keeping features with N >= {MIN_SAMPLES} samples...")
    
    # Count the number of non-missing values (samples) for each feature (row)
    sample_counts = population_matrix_raw.notna().sum(axis=1)
    
    # Filter the matrix
    population_matrix_filtered = population_matrix_raw[sample_counts >= MIN_SAMPLES]
    
    print(f"Filtered matrix shape (Features x Individuals): {population_matrix_filtered.shape}")

    # 3. Population-Level Feature Selection (Li Strategy: Highest Median Raw ER)
    
    print("Performing feature selection: Highest median raw ER per (Gene, CellType)...")
    
    # Calculate the median editing ratio across ALL individuals (axis=1) for each feature
    population_matrix_filtered['Median_Raw_ER_Population'] = population_matrix_filtered.median(axis=1)

    # Reset index to access Gene and CellType columns for grouping
    df_selection = population_matrix_filtered.reset_index()

    # Find the row index (idx) that has the maximum median within each (Gene, CellType) group
    idx_max = df_selection.groupby(['Phase3_Gene', 'CellType'])['Median_Raw_ER_Population'].idxmax()

    # Select the representative rows using these indices
    final_features_df = df_selection.loc[idx_max]

    # --- 4. Final Matrix Restructuring ---
    
    # Create the final feature identifier (e.g., APP__Bcell)
    final_features_df['FeatureID'] = final_features_df['Phase3_Gene'] + '__' + final_features_df['CellType']
    
    # Drop the temporary grouping/metadata columns
    final_features_df = final_features_df.drop(columns=['SiteID', 'Phase3_Gene', 'CellType', 'Median_Raw_ER_Population'])
    
    # Set the FeatureID as the index. The columns are still Individual IDs.
    final_edQTL_matrix = final_features_df.set_index('FeatureID')
    
    # The final edQTL matrix must have Features as rows and Individuals as columns (N_Features x N_Individuals)
    final_edQTL_matrix = final_edQTL_matrix.T.T
    
    # Ensure the index is named correctly
    final_edQTL_matrix.index.name = 'FeatureID'

    # 5. Final Output
    print(f"Final edQTL Feature Matrix shape (Features x Individuals): {final_edQTL_matrix.shape}")
    print(f"Saving final matrix to: {args.output_file}")
    
    # Save the matrix. The values are the raw editing ratios, ready for external INT in Phase 6.
    final_edQTL_matrix.to_csv(args.output_file, sep='\t', index=True, na_rep='NA')

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 5: Population-Level Collation and Feature Selection for edQTL mapping (Li Strategy).")
    parser.add_argument("--input_dir", required=True, help="Directory containing all Phase 4 output matrices.")
    parser.add_argument("--file_pattern", default="*_final_editing_matrix_p4.tsv", help="File pattern to match Phase 4 matrices.")
    parser.add_argument("--output_file", required=True, help="Path to save the single, final edQTL feature matrix.")
    parser.add_argument("--min_samples", type=int, required=True, help="Minimum number of non-missing samples required for a feature to be retained (e.g., 70).")
    args = parser.parse_args()
    
    try:
        run_phase5_collation_and_selection(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 5 Collation: {e}", file=sys.stderr)
        sys.exit(1)
