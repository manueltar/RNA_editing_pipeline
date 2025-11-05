#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
import os
import sys

# --- Main Function ---
def collate_aei_results(args):
    """
    1. Reads all individual AEI output files from the REDItools array job.
    2. Extracts the 'G-A' index (A->G on the plus strand) as the AEI.
    3. Collates results into a single (Individual x CellType) matrix.
    """
    
    input_pattern = os.path.join(args.input_dir, args.file_pattern)
    input_files = glob.glob(input_pattern)
    
    if not input_files:
        print(f"FATAL ERROR: No AEI files found matching pattern: {input_pattern}", file=sys.stderr)
        return

    print(f"Found {len(input_files)} AEI files for collation.")
    aei_list = []

    # 1. Iterate and Extract AEI
    for i, file_path in enumerate(input_files):
        try:
            # Filename format: IID_CT_ALU_ONLY.aei.tsv
            filename = os.path.basename(file_path)
            
            # Extract Individual ID (IID) and Cell Type (CT)
            # Filenames use the format: IID_CT_ALU_ONLY.aei.tsv
            parts = filename.split('_')
            individual_id = parts[0]
            cell_type = parts[1]

            # Read the AEI output file (two columns: Substitution, Index)
            # The reditools index output is two columns, tab-separated, no header
            df_aei = pd.read_csv(file_path, sep='\t', header=None, names=['Substitution', 'Index'])
            
            # Extract the A->G index, which is typically labeled 'G-A' in REDItools
            # This is the canonical AEI value.
            aei_row = df_aei[df_aei['Substitution'] == 'G-A']
            
            if aei_row.empty:
                # Should not happen if data is correctly generated, but handles edge case
                print(f"WARNING: 'G-A' substitution not found in {filename}. Skipping.", file=sys.stderr)
                continue
                
            aei_value = aei_row['Index'].iloc[0]

            aei_list.append({
                'Individual_ID': individual_id,
                'Cell_Type': cell_type,
                'AEI': aei_value
            })

        except Exception as e:
            print(f"WARNING: Skipping file {filename} due to error: {e}", file=sys.stderr)

    if not aei_list:
        print("No valid AEI data processed. Exiting.", file=sys.stderr)
        return

    # 2. Collate and Pivot to Covariate Format
    df_raw = pd.DataFrame(aei_list)
    
    # Create the final covariate column name: AEI_CellType
    df_raw['Covariate_Name'] = 'AEI_' + df_raw['Cell_Type']

    # Pivot the table: Individual_ID as Index, Covariate_Name as Columns, AEI as Values
    df_aei_final = df_raw.pivot(
        index='Individual_ID', 
        columns='Covariate_Name', 
        values='AEI'
    )
    
    df_aei_final.index.name = 'Individual_ID'

    # 3. Final Output
    print(f"Final AEI Covariate Matrix shape: {df_aei_final.shape}")
    print(f"Saving AEI covariate matrix to: {args.output_file}")
    
    # Save the matrix, ready to be merged with other covariates in Phase 6
    # NA values (if any) are written as 'NA'
    df_aei_final.to_csv(args.output_file, sep='\t', index=True, na_rep='NA')

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AEI_calculation: Collate Alu Editing Index (AEI) results.")
    parser.add_argument("--input_dir", required=True, help="Directory containing all individual AEI output files.")
    parser.add_argument("--file_pattern", default="*.aei.tsv", help="File pattern to match AEI files.")
    parser.add_argument("--output_file", required=True, help="Path to save the single, final AEI covariate matrix.")
    args = parser.parse_args()
    
    try:
        collate_aei_results(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in AEI Collation: {e}", file=sys.stderr)
        sys.exit(1)
