#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 5: Final Aggregation and Gene-Level Summarization.")
parser.add_argument("--input_matrix", required=True, help="Input matrix from Phase 4 (final_editing_matrix_p4.tsv).")
parser.add_argument("--output_file", required=True, help="Path to save the final gene-by-cell summary matrix.")
parser.add_argument("--group_by", default="gene_name", help="Column to group by for aggregation (e.g., gene_name, transcript_id).")
args = parser.parse_args()


def run_phase5_aggregation(args):
    """Loads the Phase 4 matrix, filters, aggregates by gene, and saves the summary."""
    
    print("Loading Phase 4 quantification matrix...")
    try:
        df = pd.read_csv(args.input_matrix, sep='\t')
    except Exception as e:
        print(f"FATAL ERROR: Could not load input matrix {args.input_matrix}: {e}", file=sys.stderr)
        sys.exit(1)
        
    # 1. Filter: Retain only sites that passed all global filters (VCF, SJ)
    # Sites that failed VCF/SJ have a GlobalFilterStatus != 'PASS'
    df_filtered = df[df['GlobalFilterStatus'] == 'PASS'].copy()
    print(f"Filtered matrix to {len(df_filtered)} rows (sites passed VCF/SJ filters).")

    # 2. Identify the columns to aggregate (Editing Ratios only)
    er_cols = [col for col in df_filtered.columns if col.endswith('_ER')]
    
    if not er_cols:
        print("FATAL ERROR: Could not find Editing Ratio (_ER) columns. Check column naming.", file=sys.stderr)
        sys.exit(1)

    # 3. Select relevant columns: Grouping key and ER columns
    if args.group_by not in df_filtered.columns:
        print(f"FATAL ERROR: Grouping column '{args.group_by}' not found in the input matrix.", file=sys.stderr)
        sys.exit(1)
        
    cols_to_keep = [args.group_by] + er_cols
    df_summary = df_filtered.loc[:, cols_to_keep]

    # Convert 'NA' (masked low-coverage sites from Phase 4) to true NaN for correct mean calculation
    for col in er_cols:
        # Use errors='coerce' to turn 'NA' strings into NaN. NaN values are ignored by .mean().
        df_summary[col] = pd.to_numeric(df_summary[col], errors='coerce') 

    # 4. Aggregate: Group by the chosen column and calculate the mean of all ER columns
    print(f"Aggregating by column: {args.group_by}...")
    
    gene_summary = df_summary.groupby(args.group_by)[er_cols].mean().reset_index()

    # Rename the aggregated columns for clarity (e.g., 'scRNA_file_XXX_ER' -> 'scRNA_file_XXX_AvgER')
    rename_map = {col: col.replace('_ER', '_AvgER') for col in er_cols}
    gene_summary.rename(columns=rename_map, inplace=True)
    
    print(f"Aggregation complete. Final matrix size: {len(gene_summary)} genes.")
    
    # 5. Final Output
    gene_summary.to_csv(args.output_file, sep='\t', index=False)
    print(f"Gene summary saved to {args.output_file}")


# --- 3. Main Execution ---
if __name__ == "__main__":
    run_phase5_aggregation(args)
