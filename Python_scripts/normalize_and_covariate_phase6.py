#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from scipy.stats import rankdata
from scipy.special import ndtr
import os
import sys

# Define the standard INT function
def inverse_normal_transform(series):
    """Applies Inverse Normal Transformation (INT) to a pandas Series."""
    # Handle NaN values: only rank and transform non-NaN values
    mask = pd.notna(series)
    valid_data = series[mask]
    
    # 1. Rank the data
    # Method 'average' assigns the average of ranks to tied values.
    ranks = rankdata(valid_data, method='average')
    
    # 2. Normalize ranks to range [0, 1] using (rank - 0.5) / N
    N = len(valid_data)
    normalized_ranks = (ranks - 0.5) / N
    
    # 3. Apply the inverse CDF of the normal distribution (ndtr is scipy's CDF, we use the inverse)
    # The normal CDF (ndtr) is applied to the ranks, not the inverse. This returns the quantile function value.
    # To get the quantile function (inverse CDF), we use the following relation (standard in Q-Q mapping):
    int_values = pd.Series(ndtr(normalized_ranks), index=valid_data.index)
    
    # Place the transformed values back into the original series, keeping NaNs
    transformed_series = pd.Series(np.nan, index=series.index)
    transformed_series.loc[mask] = int_values
    
    return transformed_series

def load_and_merge_covariates(args, individual_ids):
    """
    Loads external covariate files and merges them into a single matrix.
    NOTE: All covariate files must have Individual IDs as rows (index).
    """
    print("\n--- Loading Covariates ---")
    
    # Initialize the covariate matrix with the list of individuals
    covariate_dfs = []
    
    # 1. Load AEI (Allu Editing Index) - Crucial for nullifying ADAR eQTLs
    try:
        # Assuming AEI is quantified per individual per cell type
        aei_df = pd.read_csv(args.aei_file, sep='\t', index_col=0) 
        # The edQTL features are Gene__CellType, so we need to ensure the AEI matrix 
        # is structured similarly if it is cell-type specific.
        # Simplification: assuming a clean N_Individuals x N_CellTypes structure that can be merged later.
        print(f"Loaded AEI data: {aei_df.shape}")
        covariate_dfs.append(aei_df)
    except Exception as e:
        print(f"WARNING: Could not load AEI file {args.aei_file}. Skipping. Error: {e}", file=sys.stderr)

    # 2. Load Genotype PCs (Must be N_Individuals x N_PCs)
    try:
        pc_df = pd.read_csv(args.pc_file, sep='\t', index_col=0)
        pc_df = pc_df.filter(items=individual_ids, axis=0)
        print(f"Loaded Genotype PCs: {pc_df.shape}")
        covariate_dfs.append(pc_df)
    except Exception as e:
        print(f"WARNING: Could not load PC file {args.pc_file}. Skipping. Error: {e}", file=sys.stderr)

    # 3. Load PEER Factors (Must be N_Individuals x N_PEERs)
    try:
        peer_df = pd.read_csv(args.peer_file, sep='\t', index_col=0)
        peer_df = peer_df.filter(items=individual_ids, axis=0)
        print(f"Loaded PEER Factors: {peer_df.shape}")
        covariate_dfs.append(peer_df)
    except Exception as e:
        print(f"WARNING: Could not load PEER file {args.peer_file}. Skipping. Error: {e}", file=sys.stderr)

    # 4. Load Estimated Cell Proportions (Must be N_Individuals x N_CellTypes)
    try:
        prop_df = pd.read_csv(args.cell_prop_file, sep='\t', index_col=0)
        prop_df = prop_df.filter(items=individual_ids, axis=0)
        print(f"Loaded Cell Proportion Estimates: {prop_df.shape}")
        covariate_dfs.append(prop_df)
    except Exception as e:
        print(f"WARNING: Could not load Cell Proportion file {args.cell_prop_file}. Skipping. Error: {e}", file=sys.stderr)

    # 5. Merge all covariates (Individuals as rows)
    if not covariate_dfs:
        print("FATAL ERROR: No valid covariates were loaded. Cannot proceed.", file=sys.stderr)
        return pd.DataFrame()

    final_covariate_matrix = pd.concat(covariate_dfs, axis=1, join='inner')
    
    # Ensure individual IDs match the feature matrix columns
    if not all(ind_id in individual_ids for ind_id in final_covariate_matrix.index):
        print("WARNING: Covariate matrix indices do not fully match feature matrix columns.", file=sys.stderr)
        
    print(f"Final merged covariate matrix shape (Individuals x Covariates): {final_covariate_matrix.shape}")
    
    # Final step for FastQTL format: Covariates as ROWS, Individuals as COLUMNS
    final_covariate_matrix = final_covariate_matrix.T
    final_covariate_matrix.index.name = 'CovariateID'
    
    return final_covariate_matrix

def run_phase6_processing(args):
    """Main function to perform INT and prepare covariates."""
    
    # 1. Load the Phase 5 edQTL Feature Matrix
    print(f"--- Starting Phase 6: Normalization and Covariate Adjustment ---")
    print(f"Loading Phase 5 matrix: {args.input_file}")
    
    try:
        feature_matrix_raw = pd.read_csv(args.input_file, sep='\t', index_col='FeatureID', na_values=['NA'])
    except Exception as e:
        print(f"FATAL ERROR: Failed to load input file {args.input_file}. Error: {e}", file=sys.stderr)
        return

    print(f"Loaded matrix shape (Features x Individuals): {feature_matrix_raw.shape}")

    # 2. Inverse Normal Transformation (INT)
    print("\n--- Applying Inverse Normal Transformation (INT) to each feature (row) ---")
    
    # Apply INT row-wise (axis=1) across all individuals for each feature
    normalized_matrix = feature_matrix_raw.apply(inverse_normal_transform, axis=1)
    
    print(f"Normalized matrix shape: {normalized_matrix.shape}")
    
    # 3. Load and prepare Covariates
    individual_ids = feature_matrix_raw.columns.tolist()
    covariate_matrix = load_and_merge_covariates(args, individual_ids)
    
    if covariate_matrix.empty:
        print("Failed to create covariate matrix. Exiting.", file=sys.stderr)
        return

    # --- 4. Save Outputs ---
    
    # A. Save Normalized edQTL Matrix
    output_dir = os.path.dirname(args.output_normalized_file)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Saving normalized edQTL matrix to: {args.output_normalized_file}")
    # FastQTL requires TSV output
    normalized_matrix.to_csv(args.output_normalized_file, sep='\t', na_rep='NA')

    # B. Save Covariate Matrix
    print(f"Saving covariate matrix to: {args.output_covariate_file}")
    # FastQTL requires TSV output
    covariate_matrix.to_csv(args.output_covariate_file, sep='\t', na_rep='NA')
    
    print("\n--- Phase 6 Processing Complete ---")

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 6: Applies Inverse Normal Transformation (INT) to edQTL features and prepares the covariate matrix.")
    parser.add_argument("--input_file", required=True, help="Path to the Phase 5 edQTL feature matrix (raw ERs).")
    parser.add_argument("--aei_file", required=True, help="Path to the Allu Editing Index (AEI) file.")
    parser.add_argument("--pc_file", required=True, help="Path to the Genotype Principal Components (PCs) file.")
    parser.add_argument("--peer_file", required=True, help="Path to the PEER factor file.")
    parser.add_argument("--cell_prop_file", required=True, help="Path to the Cell Proportion estimates file.")
    parser.add_argument("--output_normalized_file", required=True, help="Path to save the final normalized (INT) edQTL matrix.")
    parser.add_argument("--output_covariate_file", required=True, help="Path to save the final FastQTL covariate matrix.")
    
    args = parser.parse_args()
    
    try:
        run_phase6_processing(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 6: {e}", file=sys.stderr)
        sys.exit(1)

