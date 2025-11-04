#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
import numpy as np
from scipy.stats import norm

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 6: Trait Normalization and Filtering (Inverse Normal Transformation).")
parser.add_argument("--input_matrix", required=True, help="Input gene-level AvgER matrix from Phase 5.")
parser.add_argument("--output_file", required=True, help="Path to save the final normalized trait matrix.")
parser.add_argument("--min_quantifiable_samples", type=int, default=15, help="Minimum number of samples (cell types) a gene must be quantified in to be retained.")
args = parser.parse_args()


# --- Core Function: Inverse Normal Transformation ---
def inverse_normal_transformation(series):
    """
    Applies the rank-based Inverse Normal Transformation (INT) to a pandas Series.
    The INT converts the trait values into Z-scores corresponding to the standard 
    normal distribution, ensuring the trait is normally distributed across samples.
    """
    # 1. Handle NaNs (only count non-NA values for ranking)
    non_na_values = series.dropna()
    N = len(non_na_values)
    
    if N < 2:
        # Cannot transform if too few non-NA samples
        return pd.Series(np.nan, index=series.index)

    # 2. Rank the non-NA values
    # method='average' handles ties; na_option='keep' ensures NaNs stay NaN
    ranked = series.rank(method='average', na_option='keep')
    
    # 3. Calculate the rank-based percentile (Fractional Rank)
    percentiles = (ranked - 0.5) / N
    
    # 4. Apply the Inverse Normal (Quantile) function
    transformed = pd.Series(norm.ppf(percentiles.dropna()), index=non_na_values.index)
    
    # 5. Combine with original series index to reintroduce NaNs where needed
    result = pd.Series(np.nan, index=series.index)
    result.update(transformed)
    
    return result


def run_phase6_normalization(args):
    """Loads the Phase 5 matrix, filters for sample coverage, applies INT, and saves the final trait matrix."""

    print("Loading Phase 5 gene summary matrix...")
    try:
        df = pd.read_csv(args.input_matrix, sep='\t')
    except Exception as e:
        print(f"FATAL ERROR: Could not load input matrix {args.input_matrix}: {e}", file=sys.stderr)
        sys.exit(1)
        
    # --- A. Preparation ---
    # Gene names are in the first column, treat the rest as trait columns
    gene_col = df.columns[0]
    trait_cols = df.columns[1:]
    
    # Convert 'NA' strings (if any survived Phase 5) to true NaN, and ensure floats
    df[trait_cols] = df[trait_cols].apply(pd.to_numeric, errors='coerce')

    # --- B. Trait Filtering (Required for EdQTL Power) ---
    
    # Count the number of non-NaN values (quantified samples) for each gene (row)
    df['Quantified_Samples_Count'] = df[trait_cols].count(axis=1) 
    
    df_filtered = df[df['Quantified_Samples_Count'] >= args.min_quantifiable_samples].copy()
    
    print(f"Filtered matrix from {len(df)} to {len(df_filtered)} genes (retaining those with >= {args.min_quantifiable_samples} quantified samples).")
    
    if df_filtered.empty:
        print("FATAL ERROR: No genes passed the minimum sample quantification filter.", file=sys.stderr)
        sys.exit(1)

    # --- C. Inverse Normal Transformation (INT) ---
    # Transpose the data to make samples rows and genes columns (easier to apply INT across columns/traits)
    df_transposed = df_filtered.set_index(gene_col)[trait_cols].T
    
    print(f"Applying Inverse Normal Transformation to {len(df_transposed.columns)} traits...")

    # Apply INT to each column (which represents one gene's AvgER values across samples)
    df_normalized_transposed = df_transposed.apply(inverse_normal_transformation, axis=0)

    # Transpose back to original Gene-by-Cell format
    df_normalized = df_normalized_transposed.T.reset_index()
    
    # Rename columns to reflect normalization
    new_cols = {col: col.replace('_AvgER', '_INT') for col in df_normalized.columns if col.endswith('_AvgER')}
    df_normalized.rename(columns=new_cols, inplace=True)
    
    # --- D. Final Output ---
    df_normalized.to_csv(args.output_file, sep='\t', index=False)
    print(f"Final normalized trait matrix saved to {args.output_file}")


# --- 3. Main Execution ---
if __name__ == "__main__":
    run_phase6_normalization(args)
