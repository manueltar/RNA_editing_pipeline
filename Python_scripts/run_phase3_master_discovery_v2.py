#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob
import sys

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 3 (Reworked): Master Site Discovery and Consensus Filtering.")
parser.add_argument("--raw_calls_dir", required=True, help="Directory containing all raw RED-ML and REDItools chunk outputs.")
parser.add_argument("--output_file", required=True, help="Path to save the final consensus master site list (Chr, Pos, Ref, Alt only).")
parser.add_argument("--min_edit_level", type=float, default=0.1, help="Minimum Editing Level (e.g., 0.1 = 10%%) required for a site to be kept by a tool.")
args = parser.parse_args()

# --- 2. Site Discovery and Aggregation Logic ---

def load_and_aggregate_raw_calls(raw_calls_dir, min_edit_level):
    """Loads and aggregates all raw call files, applying the minimum Edit Level filter."""

    all_sites = []
    
    redml_files = glob.glob(os.path.join(raw_calls_dir, 'raw_redml_output_chunk_*.tsv'))
    reditools_files = glob.glob(os.path.join(raw_calls_dir, 'raw_reditools_output_chunk_*.tsv'))

    if not redml_files or not reditools_files:
        print("ERROR: Could not find raw call files. Check raw_calls_dir path and file patterns.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(redml_files)} RED-ML files and {len(reditools_files)} REDItools files.")
    
    # --- A. Process RED-ML Files (Calculate and filter by Edit Level) ---
    for f in redml_files:
        try:
            # We assume RED-ML output has 'VariantReads' and 'TotalReads'
            df = pd.read_csv(f, sep='\t', usecols=['Chr', 'Pos', 'Ref', 'Alt', 'VariantReads', 'TotalReads'])
            
            # Calculate Edit Level (ER)
            df['EditLevel'] = df['VariantReads'] / df['TotalReads'].replace(0, pd.NA) 
            df = df.dropna(subset=['EditLevel']) 

            # **EXPLICIT QUANTITATIVE QC FILTER**
            df = df[df['EditLevel'] >= min_edit_level].copy()
            
            df['Tool'] = 'REDML'
            all_sites.append(df.loc[:, ['Chr', 'Pos', 'Ref', 'Alt', 'Tool']])
        except Exception as e:
            print(f"Warning: Failed to process RED-ML file {os.path.basename(f)}: {e}", file=sys.stderr)

    # --- B. Process REDItools Files (Filter by Edit Level) ---
    col_map = {'Chromosome': 'Chr', 'Position': 'Pos', 'Reference': 'Ref', 'Edited': 'Alt', 'Editing_Level': 'EditLevel'}

    for f in reditools_files:
        try:
            # Assuming REDItools provides an 'Editing_Level' column
            df = pd.read_csv(f, sep='\t').rename(columns=col_map)
            df = df.loc[:, ['Chr', 'Pos', 'Ref', 'Alt', 'EditLevel']] 
            df = df.dropna(subset=['EditLevel'])

            # **EXPLICIT QUANTITATIVE QC FILTER**
            df = df[df['EditLevel'] >= min_edit_level].copy() 

            df['Tool'] = 'REDItools'
            all_sites.append(df.loc[:, ['Chr', 'Pos', 'Ref', 'Alt', 'Tool']])
        except Exception as e:
            print(f"Warning: Failed to process REDItools file {os.path.basename(f)}: {e}", file=sys.stderr)

    if not all_sites:
        print("FATAL ERROR: No valid site data was loaded after filtering.", file=sys.stderr)
        sys.exit(1)

    return pd.concat(all_sites, ignore_index=True)


def run_master_discovery(args):
    """Runs the full discovery and minimal filtering process."""

    raw_df = load_and_aggregate_raw_calls(args.raw_calls_dir, args.min_edit_level)
    initial_count = len(raw_df)
    print(f"\nTotal raw entries loaded after Edit Level >= {args.min_edit_level} filter: {initial_count}")

    raw_df['SiteID'] = raw_df['Chr'].astype(str) + ':' + raw_df['Pos'].astype(str) + ':' + raw_df['Ref'] + '>' + raw_df['Alt']

    # --- C. Canonical Editing Filter (A>G or T>C) ---
    canonical_df = raw_df[
        ((raw_df['Ref'] == 'A') & (raw_df['Alt'] == 'G')) |
        ((raw_df['Ref'] == 'T') & (raw_df['Alt'] == 'C'))
    ].copy()

    print(f"Sites after Canonical Editing Filter (A>G or T>C): {len(canonical_df)}")

    # --- D. Consensus Filter (Found by both tools) ---
    consensus_check = canonical_df.groupby('SiteID')['Tool'].nunique().reset_index(name='ToolCount')
    consensus_sites = consensus_check[consensus_check['ToolCount'] == 2]['SiteID']

    final_master_df = canonical_df[canonical_df['SiteID'].isin(consensus_sites)].copy()
    final_master_df = final_master_df[['Chr', 'Pos', 'Ref', 'Alt']].drop_duplicates()

    print(f"Sites after Consensus Filter (Found by both RED-ML & REDItools): {len(final_master_df)}")

    # --- E. Final Output ---
    final_master_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"\nPhase 3 Master Site Discovery complete. Unique sites written to: {args.output_file}")


# --- 3. Main Execution ---
if __name__ == "__main__":
    try:
        run_master_discovery(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 3 Master Discovery: {e}", file=sys.stderr)
        sys.exit(1)
