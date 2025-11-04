#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob
import sys
from io import StringIO
from itertools import chain

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 3: Master Site Discovery and Consensus Filtering.")
# CRITICAL FIX: Match the argument name from the shell script
parser.add_argument("--root_search_dir", required=True, 
                    help="Top-level directory to recursively search for raw call files (IND_XXXX/P2_.../).")
parser.add_argument("--output_file", required=True, 
                    help="Path to save the final consensus master site list (Chr, Pos, Ref, Alt only).")
parser.add_argument("--min_edit_level", type=float, default=0.1, 
                    help="Minimum Editing Level (e.1. 0.1 = 10%%) required for a site to be kept by REDItools.")
parser.add_argument("--threads", type=int, default=4,
                    help="Number of threads for parallel processing (used for I/O optimization).") # Added thread argument
args = parser.parse_args()

# --- 2. Helper Function for REDItools Substitution Parsing ---

def explode_reditools_substitutions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Takes a REDItools DataFrame and normalizes the 'AllSubs' column 
    by splitting space-separated substitutions (e.g., 'A>G C>T') into 
    separate rows, applying the Editing Level to each.
    """
    # 1. Map Ref and AllSubs into a list of (Ref, Alt, EditLevel) for each row
    sub_data = []
    for _, row in df.iterrows():
        # AllSubs is like "A>G C>T", Frequency is a float
        if row['AllSubs'] and row['AllSubs'] != '-':
            substitutions = row['AllSubs'].split(' ')
            for sub in substitutions:
                # The assumption here is that if a position has multiple substitutions, 
                # the single 'Frequency' (EditLevel) reported applies to the sum, 
                # or the tool reports the highest frequency substitution.
                # For safety, we keep the original frequency but derive the Alt base.
                if len(sub) == 3 and sub[1] == '>':
                    sub_data.append({
                        'Chr': row['Chr'],
                        'Pos': row['Pos'],
                        'Ref': row['Reference'],
                        'Alt': sub[2],
                        'EditLevel': row['Frequency']
                    })
    
    if not sub_data:
        return pd.DataFrame()
        
    return pd.DataFrame(sub_data)


# --- 3. Site Discovery and Aggregation Logic ---

def load_and_aggregate_raw_calls(root_search_dir, min_edit_level):
    """Loads and aggregates all raw call files with recursive search."""

    # CRITICAL FIX: Use recursive globbing and correct final file names
    redml_files = glob.glob(os.path.join(root_search_dir, '**', '*_redml_raw.tsv'), recursive=True)
    reditools_files = glob.glob(os.path.join(root_search_dir, '**', '*_reditools_raw.tsv'), recursive=True)

    if not redml_files or not reditools_files:
        print("ERROR: Could not find raw call files. Check root_search_dir path and file patterns.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(redml_files)} RED-ML files and {len(reditools_files)} REDItools files.")
    
    all_sites = []
    
    # --- A. Process RED-ML Files (Calculate and filter by Edit Level) ---
    print("Processing RED-ML files...")
    redml_columns = ['Chr', 'Pos', 'Ref', 'Alt', 'VariantReads', 'TotalReads']
    for f in redml_files:
        try:
            # We assume the RED-ML output is now standardized
            df = pd.read_csv(f, sep='\t', usecols=redml_columns, dtype={'Chr': str, 'Pos': int})
            
            # Calculate Edit Level (ER). Note: This is a filter for RED-ML as well.
            # RED-ML output is assumed to be filtered for sites, but we apply the ER filter
            # for consistency with the REDItools filtering step.
            df['EditLevel'] = df['VariantReads'] / df['TotalReads'].replace(0, pd.NA) 
            df = df.dropna(subset=['EditLevel']) 

            # **EXPLICIT QUANTITATIVE QC FILTER**
            df = df[df['EditLevel'] >= min_edit_level].copy()
            
            df['Tool'] = 'REDML'
            # Extract sample ID from filename for possible future per-sample consensus
            df['SampleID'] = os.path.basename(f).split('_redml_raw.tsv')[0]
            all_sites.append(df.loc[:, ['Chr', 'Pos', 'Ref', 'Alt', 'Tool', 'SampleID']])
        except Exception as e:
            # If a file is malformed, we log it and continue
            print(f"Warning: Failed to process RED-ML file {os.path.basename(f)}: {e}", file=sys.stderr)

    # --- B. Process REDItools Files (Filter by Edit Level and Explode) ---
    print("Processing REDItools files...")
    reditools_columns = ['Region', 'Position', 'Reference', 'AllSubs', 'Frequency']
    
    for f in reditools_files:
        try:
            # CRITICAL FIX: Use the correct REDItools3 column names
            df = pd.read_csv(f, sep='\t', usecols=reditools_columns, dtype={'Region': str, 'Position': int})
            
            # Map REDItools columns to standardized names
            df.rename(columns={'Region': 'Chr', 'Position': 'Pos'}, inplace=True)
            
            # CRITICAL FIX: Explode and normalize the 'AllSubs' column to get 'Alt'
            df_exploded = explode_reditools_substitutions(df)
            
            # **EXPLICIT QUANTITATIVE QC FILTER**
            df_exploded = df_exploded[df_exploded['EditLevel'] >= min_edit_level].copy() 

            df_exploded['Tool'] = 'REDItools'
            # Extract sample ID from filename for possible future per-sample consensus
            df_exploded['SampleID'] = os.path.basename(f).split('_reditools_raw.tsv')[0]
            all_sites.append(df_exploded.loc[:, ['Chr', 'Pos', 'Ref', 'Alt', 'Tool', 'SampleID']])
            
        except Exception as e:
            print(f"Warning: Failed to process REDItools file {os.path.basename(f)}: {e}", file=sys.stderr)

    if not all_sites:
        print("FATAL ERROR: No valid site data was loaded after filtering.", file=sys.stderr)
        sys.exit(1)

    return pd.concat(all_sites, ignore_index=True)


def run_master_discovery(args):
    """Runs the full discovery and minimal filtering process."""

    # CRITICAL FIX: Pass the correct argument name
    raw_df = load_and_aggregate_raw_calls(args.root_search_dir, args.min_edit_level)
    initial_count = len(raw_df)
    print(f"\nTotal raw entries loaded after Edit Level >= {args.min_edit_level} filter: {initial_count}")

    # Use a robust key that uniquely identifies the site
    raw_df['SiteID'] = raw_df['Chr'].astype(str) + ':' + raw_df['Pos'].astype(str) + ':' + raw_df['Ref'] + '>' + raw_df['Alt']

    # --- C. Canonical Editing Filter (A>G or T>C) ---
    canonical_df = raw_df[
        ((raw_df['Ref'] == 'A') & (raw_df['Alt'] == 'G')) |
        ((raw_df['Ref'] == 'T') & (raw_df['Alt'] == 'C'))
    ].copy()

    print(f"Sites after Canonical Editing Filter (A>G or T>C): {len(canonical_df)}")

    # --- D. Consensus Filter (Found by both tools somewhere in the project) ---
    # This checks for project-wide consensus, not per-sample consensus.
    consensus_check = canonical_df.groupby('SiteID')['Tool'].nunique().reset_index(name='ToolCount')
    consensus_sites = consensus_check[consensus_check['ToolCount'] == 2]['SiteID']

    final_master_df = canonical_df[canonical_df['SiteID'].isin(consensus_sites)].copy()
    final_master_df = final_master_df[['Chr', 'Pos', 'Ref', 'Alt']].drop_duplicates()

    print(f"Sites after Project-Wide Consensus Filter (Found by both RED-ML & REDItools): {len(final_master_df)}")

    # --- E. Final Output ---
    final_master_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"\nPhase 3 Master Site Discovery complete. Unique sites written to: {args.output_file}")


# --- 4. Main Execution ---
if __name__ == "__main__":
    try:
        run_master_discovery(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 3 Master Discovery: {e}", file=sys.stderr)
        sys.exit(1)
