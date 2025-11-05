#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from statsmodels.sandbox.stats.multicomp import multipletests
import sys
import os

# --- Main Function ---
def process_fastqtl_results(args):
    """
    1. Loads FastQTL output.
    2. Identifies the best permutation P-value (P_beta) for each feature.
    3. Calculates Q-values (FDR correction) based on P_beta.
    4. Filters and saves the final list of significant lead edQTLs.
    """
    FDR_THRESHOLD = args.fdr_threshold # Default: 0.05

    print(f"--- Starting Phase 8: Processing FastQTL Results ---")
    print(f"Loading FastQTL results from: {args.input_file}")
    
    # 1. Load Data
    try:
        # FastQTL output is space-separated, compressed
        df_raw = pd.read_csv(args.input_file, sep='\s+', compression='infer', header=None)
    except FileNotFoundError:
        print(f"FATAL ERROR: Input file not found at {args.input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"FATAL ERROR while loading file: {e}", file=sys.stderr)
        sys.exit(1)

    # FastQTL Output Columns (assuming chunk 1 1 output):
    # 0: FeatureID
    # 1: n_snps
    # 2: P_nominal (best nominal p-value)
    # 3: SNP_ID (of best nominal)
    # 4: SNP_distance
    # 5: nominal_P_value_beta (FastQTL-specific: used for P_perm calculation)
    # 6: P_permutation (empirical p-value, P_perm)
    # 7: P_beta (Best empirical p-value, P_beta)
    
    # Use the column names that correspond to the permutation-based output
    df_raw.columns = ['FeatureID', 'n_snps', 'P_nominal', 'SNP_ID', 'SNP_distance', 
                      'nominal_P_value_beta', 'P_permutation', 'P_beta']
    
    # Check for empty results
    if df_raw.empty:
        print("WARNING: FastQTL output file is empty. Exiting.", file=sys.stderr)
        sys.exit(0)

    # 2. FDR Correction (Q-value Calculation)
    print(f"Total features tested: {len(df_raw)}")
    
    # The P_beta column (Best empirical P-value after beta-approximation) is used for FDR correction
    p_values = df_raw['P_beta'].values
    
    # Apply Benjamini-Hochberg procedure
    reject, q_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Add Q-values and Significance status to the dataframe
    df_raw['Q_value'] = q_values
    df_raw['Significant'] = reject
    
    print(f"Number of features passing FDR < {FDR_THRESHOLD} (Q-value): {np.sum(reject)}")

    # 3. Filtering and Final Output
    df_final_edQTLs = df_raw[df_raw['Q_value'] < FDR_THRESHOLD].copy()

    # Select and rename columns for clarity
    df_final_edQTLs = df_final_edQTLs[['FeatureID', 'SNP_ID', 'P_nominal', 'P_beta', 'Q_value', 'SNP_distance']]
    df_final_edQTLs.rename(columns={
        'SNP_ID': 'Lead_SNP_ID',
        'P_nominal': 'P_Nominal_Best',
        'P_beta': 'P_Empirical_Best',
        'Q_value': 'FDR_Q_value',
        'SNP_distance': 'Lead_SNP_Distance'
    }, inplace=True)
    
    # Save the full table and the significant subset
    
    # Save full results with Q-values
    full_output_file = args.output_file.replace('.tsv', '_full_with_qval.tsv')
    df_raw.to_csv(full_output_file, sep='\t', index=False)
    print(f"Saved full results (with Q-values) to: {full_output_file}")
    
    # Save significant results
    print(f"Saving {len(df_final_edQTLs)} significant edQTLs to: {args.output_file}")
    df_final_edQTLs.to_csv(args.output_file, sep='\t', index=False)

    print("--- Phase 8 Processing Complete ---")

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 8: Multiple Testing Correction and Lead SNP Identification.")
    parser.add_argument("--input_file", required=True, help="Input FastQTL results file (e.g., from Phase 7).")
    parser.add_argument("--output_file", required=True, help="Path to save the final table of significant edQTLs (FDR < 0.05).")
    parser.add_argument("--fdr_threshold", type=float, default=0.05, help="FDR threshold for significance (Q-value).")
    args = parser.parse_args()
    
    try:
        # Handle compressed files (.gz) for input
        if args.input_file.endswith('.gz'):
            args.input_file = args.input_file.replace('.gz', '')
        
        process_fastqtl_results(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 8 Processing: {e}", file=sys.stderr)
        sys.exit(1)
