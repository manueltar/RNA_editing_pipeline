#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
import os
import sys
from statsmodels.sandbox.stats.multicomp import multipletests

# --- Main Function ---
def run_fdr_correction(input_files, output_dir):
    """
    Loads FastQTL results, performs FDR correction using the BH method 
    on the empirical p-values, and saves the final significant results.
    """
    
    # 1. Load Data
    all_pvalues = []
    
    # FastQTL output columns of interest:
    # 1: gene_id (Feature ID), 2: variant_id (SNP ID), 12: p_beta (Empirical P-value)
    
    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"WARNING: Input file not found: {file_path}. Skipping.", file=sys.stderr)
            continue
            
        print(f"Loading and processing: {os.path.basename(file_path)}")
        
        try:
            # Read only the necessary columns (0-based indexing)
            df = pd.read_csv(file_path, sep='\t', compression='gzip', 
                             usecols=[0, 1, 11], 
                             names=['feature_id', 'variant_id', 'p_empirical'],
                             skiprows=1) # Skip the header if present, or FastQTL summary line
            
            df['source_file'] = os.path.basename(file_path)
            all_pvalues.append(df)
        except Exception as e:
            print(f"ERROR processing {file_path}: {e}", file=sys.stderr)

    if not all_pvalues:
        print("FATAL ERROR: No valid data loaded for correction. Exiting.", file=sys.stderr)
        return
        
    df_combined = pd.concat(all_pvalues, ignore_index=True)
    
    # 2. FDR Correction (Benjamini/Hochberg)
    # We apply the correction to ALL P-values combined (main edQTLs + AEI-QTLs) 
    # as the edQTL study defines the total test burden.
    
    # Filter out potential NaNs if any issues occurred in permutation
    df_combined = df_combined.dropna(subset=['p_empirical'])
    
    # Note: FastQTL empirical p-values are typically well-calibrated, so this single BH
    # correction on the combined set is standard for a genome-wide edQTL study.
    
    print(f"Applying BH FDR correction on {len(df_combined)} total tests...")
    
    # Use the Benjamini/Hochberg procedure to control False Discovery Rate (FDR)
    reject, qvals, _, _ = multipletests(df_combined['p_empirical'], 
                                        alpha=0.05, 
                                        method='fdr_bh')
    
    df_combined['q_value'] = qvals
    df_combined['FDR_significant'] = reject

    # 3. Filtering and Lead SNP Identification
    # Identify significant results
    df_significant = df_combined[df_combined['FDR_significant']]
    print(f"Total significant edQTLs/AEI-QTLs at q < 0.05: {len(df_significant)}")
    
    # Group by Feature ID (editing site or AEI_CellType) and find the lead SNP (lowest P-value)
    # The 'idxmin' on the empirical P-value will select the lead variant
    idx = df_significant.groupby('feature_id')['p_empirical'].idxmin()
    df_lead_snps = df_significant.loc[idx]
    
    print(f"Total unique lead edQTLs/AEI-QTLs identified: {len(df_lead_snps)}")

    # 4. Save Final Output
    final_output_path = os.path.join(output_dir, "master_edQTL_AEIQTL_lead_snps.tsv")
    print(f"Saving final lead SNPs to: {final_output_path}")
    df_lead_snps.sort_values(by='q_value').to_csv(final_output_path, sep='\t', index=False)
    
    # Also save the full corrected table for later use
    full_output_path = os.path.join(output_dir, "master_edQTL_AEIQTL_full_corrected.tsv.gz")
    df_combined.to_csv(full_output_path, sep='\t', index=False, compression='gzip')
    print("Full corrected results saved.")

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 8: FDR Correction and Lead SNP Identification.")
    parser.add_argument("--input_edqtl_results", required=True, nargs='+', help="List of FastQTL edQTL result files (e.g., chr*.txt.gz).")
    parser.add_argument("--input_aeiqtl_results", required=True, nargs='+', help="List of FastQTL AEI-QTL result files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the final corrected results.")
    args = parser.parse_args()
    
    # Combine the input lists for processing
    all_inputs = args.input_edqtl_results + args.input_aeiqtl_results
    
    try:
        run_fdr_correction(all_inputs, args.output_dir)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 8: {e}", file=sys.stderr)
        sys.exit(1)
