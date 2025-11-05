#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import os

def prepare_aeiqtl_matrices(args):
    """
    1. Loads the AEI matrix (from AEI_calculation, which is Individual x AEI_CellType).
    2. Transposes it for use as the FastQTL phenotype (Feature x Individual).
    3. Loads general covariates (PCs, PEER) and saves them for AEI-QTL.
    """
    print("--- Preparing AEI-QTL Phenotype and Covariate Matrices ---")
    
    # 1. Load the AEI Matrix (which is our new phenotype)
    try:
        # AEI matrix shape: Individuals (rows) x AEI_CellType (columns)
        df_aei = pd.read_csv(args.input_aei_covariates, sep='\t', index_col='Individual_ID')
        print(f"Loaded AEI matrix (as Phenotype): {df_aei.shape}")
    except Exception as e:
        print(f"FATAL ERROR loading AEI covariate file: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Transpose for FastQTL Phenotype format
    # FastQTL requires Feature (rows) x Individual (columns)
    df_phenotype_aeiqtl = df_aei.T
    df_phenotype_aeiqtl.index.name = 'feature_id'
    
    # 3. Load Genotype PCs
    try:
        df_pcs = pd.read_csv(args.input_genotype_pcs, sep='\t', index_col='Individual_ID')
    except Exception:
        print("WARNING: Genotype PCs file not found or failed to load. Using empty matrix.", file=sys.stderr)
        df_pcs = pd.DataFrame()
        
    # 4. Load PEER Factors
    try:
        df_peer = pd.read_csv(args.input_peer_factors, sep='\t', index_col='Individual_ID')
    except Exception:
        print("WARNING: PEER Factors file not found or failed to load. Using empty matrix.", file=sys.stderr)
        df_peer = pd.DataFrame()
        
    # 5. Merge Covariates (excluding the AEI covariate itself)
    # We use a subset of individuals common to the AEI phenotype and all covariates
    individuals = df_phenotype_aeiqtl.columns.intersection(df_pcs.index).intersection(df_peer.index)
    
    df_covariates_merged = pd.concat([df_pcs, df_peer], axis=1)
    
    # Finalize phenotype and covariates on the common set of individuals
    df_phenotype_aeiqtl = df_phenotype_aeiqtl.loc[:, individuals]
    df_covariates_aeiqtl = df_covariates_merged.loc[individuals].fillna(0)
    
    df_covariates_aeiqtl.index.name = 'individual_id'
    df_covariates_aeiqtl = df_covariates_aeiqtl.loc[:, df_covariates_aeiqtl.nunique() > 1] # Remove non-variable

    # 6. Save Outputs
    print(f"Saving AEI-QTL Phenotype Matrix ({df_phenotype_aeiqtl.shape}) to: {args.output_phenotype}")
    df_phenotype_aeiqtl.to_csv(args.output_phenotype, sep='\t', index=True)

    print(f"Saving AEI-QTL Covariate Matrix ({df_covariates_aeiqtl.shape}) to: {args.output_covariates}")
    df_covariates_aeiqtl.to_csv(args.output_covariates, sep='\t', index=True)
    
    print("--- AEI-QTL Matrix Preparation Complete ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 6b: Prepare matrices for AEI-QTL mapping.")
    parser.add_argument("--input_aei_covariates", required=True, help="Input file of AEI covariates (Phase AEI_calculation output).")
    parser.add_argument("--input_genotype_pcs", required=True, help="Input file of Genotype Principal Components (PCs).")
    parser.add_argument("--input_peer_factors", required=True, help="Input file of PEER Factors.")
    parser.add_argument("--output_phenotype", required=True, help="Path to save the final AEI-QTL phenotype matrix (FastQTL input).")
    parser.add_argument("--output_covariates", required=True, help="Path to save the final AEI-QTL covariate matrix (FastQTL input).")
    args = parser.parse_args()
    prepare_aeiqtl_matrices(args)
