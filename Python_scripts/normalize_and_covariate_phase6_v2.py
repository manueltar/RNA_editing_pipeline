#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import sys
import os
from scipy.stats import rankdata
from scipy.special import ndtr

# --- Core INT Function ---
def inverse_normal_transform(data):
    """
    Applies Inverse Normal Transformation (INT) to a single array (column).
    INT converts the data to a normal distribution.
    
    Equivalent to: norm.ppf((rank(X, ties='average') - 0.5) / n)
    """
    # 1. Rank the data (average rank for ties)
    ranked_data = rankdata(data, method='average')
    
    # 2. Normalize ranks to range [0, 1]
    n = len(data)
    normalized_ranks = (ranked_data - 0.5) / n
    
    # 3. Apply the inverse of the Normal Cumulative Distribution Function (CDF)
    # ndtr is the Normal CDF, so we use it on the normalized ranks to get the transformed value
    return ndtr(normalized_ranks)

# --- Main Function ---
def run_normalization_and_merge(args):
    """
    1. Loads feature data (Phase 5 output).
    2. Loads and merges all covariate data (AEI, PCs, PEER).
    3. Applies Inverse Normal Transformation (INT) to features.
    4. Saves final phenotype and covariate matrices.
    """
    print(f"--- Starting Phase 6: Normalization and Covariate Merge ---")
    
    # 1. Load Feature Matrix (Phenotype)
    try:
        # Features are columns (Edit_Site__CellType), Individuals are rows (Index)
        df_phenotype = pd.read_csv(args.input_features, sep='\t', index_col='Individual_ID')
        print(f"Loaded feature matrix: {df_phenotype.shape}")
    except Exception as e:
        print(f"FATAL ERROR loading feature file: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Load Covariate Matrices
    covariates = []
    
    # 2a. Load AEI Covariates (from AEI_calculation phase)
    try:
        df_aei = pd.read_csv(args.input_aei_covariates, sep='\t', index_col='Individual_ID')
        df_aei = df_aei.dropna(axis=1, how='all') # Drop any columns that are all NA
        print(f"Loaded AEI covariates: {df_aei.shape}")
        covariates.append(df_aei)
    except Exception as e:
        print(f"FATAL ERROR loading AEI covariates: {e}", file=sys.stderr)
        sys.exit(1)
        
    # 2b. Load Genotype PCs (Placeholder)
    try:
        df_pcs = pd.read_csv(args.input_genotype_pcs, sep='\t', index_col='Individual_ID')
        print(f"Loaded Genotype PCs: {df_pcs.shape}")
        covariates.append(df_pcs)
    except Exception as e:
        print(f"WARNING: Genotype PCs file not found or failed to load ({e}). Skipping.")
        
    # 2c. Load PEER Factors (Placeholder - assuming K=60 factors)
    try:
        df_peer = pd.read_csv(args.input_peer_factors, sep='\t', index_col='Individual_ID')
        print(f"Loaded PEER Factors: {df_peer.shape}")
        covariates.append(df_peer)
    except Exception as e:
        print(f"WARNING: PEER Factors file not found or failed to load ({e}). Skipping.")
        
    # 3. Merge Covariates
    # Use outer merge to keep all individuals, filling missing covariate values with 0
    df_covariates_merged = pd.concat(covariates, axis=1)
    # Filter individuals to those present in the phenotype data
    individuals = df_phenotype.index.intersection(df_covariates_merged.index)
    
    df_phenotype = df_phenotype.loc[individuals]
    df_covariates_final = df_covariates_merged.loc[individuals].fillna(0) # Treat missing covariates as zero effect

    # Drop columns from the covariate matrix that have no variation (e.g., all 0)
    df_covariates_final = df_covariates_final.loc[:, df_covariates_final.nunique() > 1]
    
    print(f"Final merged and filtered Covariates: {df_covariates_final.shape}")
    print(f"Final Phenotype matrix to transform: {df_phenotype.shape}")
    
    # 4. Apply Inverse Normal Transformation (INT)
    print("Applying Inverse Normal Transformation (INT) to all feature columns...")
    
    # Apply INT column-wise (axis=0)
    df_phenotype_int = df_phenotype.apply(inverse_normal_transform, axis=0)

    # 5. Save Outputs (FastQTL format)
    
    # a. Phenotype Matrix (Features are now rows, Individuals are columns - FastQTL format)
    # The INT matrix must be transposed for FastQTL
    df_phenotype_int_fastqtl = df_phenotype_int.T
    df_phenotype_int_fastqtl.index.name = 'feature_id'
    
    print(f"Saving INT Phenotype Matrix ({df_phenotype_int_fastqtl.shape}) to: {args.output_phenotype}")
    df_phenotype_int_fastqtl.to_csv(args.output_phenotype, sep='\t', index=True)

    # b. Covariate Matrix (Individuals are rows, Covariates are columns - FastQTL format)
    df_covariates_final.index.name = 'individual_id' # Set column name for FastQTL
    
    print(f"Saving Covariate Matrix ({df_covariates_final.shape}) to: {args.output_covariates}")
    df_covariates_final.to_csv(args.output_covariates, sep='\t', index=True)

    print("--- Phase 6 Normalization and Merge Complete ---")

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 6: Normalization (INT) and Covariate Merge.")
    parser.add_argument("--input_features", required=True, help="Input file of selected features (Phase 5 output).")
    parser.add_argument("--input_aei_covariates", required=True, help="Input file of AEI covariates (AEI_calculation output).")
    parser.add_argument("--input_genotype_pcs", required=True, help="Input file of Genotype Principal Components (PCs).")
    parser.add_argument("--input_peer_factors", required=True, help="Input file of PEER Factors.")
    parser.add_argument("--output_phenotype", required=True, help="Path to save the final INT-transformed phenotype matrix (FastQTL input).")
    parser.add_argument("--output_covariates", required=True, help="Path to save the final merged covariate matrix (FastQTL input).")
    args = parser.parse_args()
    
    try:
        run_normalization_and_merge(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 6: {e}", file=sys.stderr)
        sys.exit(1)
