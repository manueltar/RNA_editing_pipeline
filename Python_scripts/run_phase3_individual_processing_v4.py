#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os
import glob
import sys
import scipy.stats 


# --- 1. Argument Parsing and Setup ---
parser = argparse.ArgumentParser(description="Phase 3: Individual Aggregation, Annotation, and Normalization.")
parser.add_argument("--individual_id", required=True, 
                    help="The unique ID of the individual being processed (e.g., IND_0001).")
parser.add_argument("--root_search_dir", required=True, 
                    help="Top-level directory to recursively search for raw call files.")
parser.add_argument("--output_file", required=True, 
                    help="Path to save the final normalized and annotated matrix for this individual.")
parser.add_argument("--min_edit_level", type=float, default=0.1, 
                    help="Minimum Editing Level (e.g., 0.1 = 10%%) required for a site to be kept.")
# CRITICAL FIX: Changed argument name to rediportal_bed
parser.add_argument("--rediportal_bed", required=True,
                    help="Path to the REDIPortal Known Sites BED file.")
parser.add_argument("--gene_map_tsv", required=True,
                    help="Path to the pre-processed Site-to-Gene map TSV file.")
args = parser.parse_args()


# --- 2. Helper Functions (Annotation & Normalization) ---

def explode_reditools_substitutions(df: pd.DataFrame) -> pd.DataFrame:
    """Explodes the REDItools 'AllSubs' column into multiple rows, standardizing columns."""
    sub_data = []
    if df.empty: return pd.DataFrame()
        
    for _, row in df.iterrows():
        if pd.notna(row['AllSubs']) and row['AllSubs'] != '-':
            substitutions = row['AllSubs'].split(' ')
            for sub in substitutions:
                if len(sub) == 3 and sub[1] == '>':
                    sub_data.append({
                        'Chr': row['Chr'],
                        'Pos': row['Pos'],
                        'Reference': row['Reference'],
                        'Alt': sub[2],
                        'EditLevel': row['Frequency']
                    })
    if not sub_data: return pd.DataFrame()
    return pd.DataFrame(sub_data).rename(columns={'Reference': 'Ref'})

def annotate_rediportal_status(df: pd.DataFrame, rediportal_bed_path: str) -> pd.DataFrame:
    """Adds a 'REDIPortal_Status' column ('Known'/'Novel')."""
    # CRITICAL FIX: Updated function and column names
    df['REDIPortal_Status'] = 'Novel'
    try:
        # Load BED file (Chr, Start, End are standard)
        rediportal_sites = pd.read_csv(rediportal_bed_path, sep='\t', header=None, usecols=[0,1], compression='infer', dtype={0: str})
        rediportal_sites.columns = ['Chr', 'Start']
        rediportal_sites['Pos'] = rediportal_sites['Start'] + 1 
        rediportal_sites['KnownID'] = rediportal_sites['Chr'].astype(str) + ':' + rediportal_sites['Pos'].astype(str)
        known_site_ids = set(rediportal_sites['KnownID'])
        
        # Match using the Chr:Pos part of the SiteID index
        df.loc[df.index.to_series().apply(lambda x: x.split(':')[0] + ':' + x.split(':')[1]).isin(known_site_ids), 'REDIPortal_Status'] = 'Known'
    except Exception as e:
        print(f"  Warning: Failed to load REDIPortal: {e}", file=sys.stderr)
    return df

def annotate_functional_region(df: pd.DataFrame, gene_map_tsv_path: str) -> pd.DataFrame:
    """Adds 'Functional_Region' and 'Gene' columns based on the Gene Map."""
    try:
        gene_map = pd.read_csv(gene_map_tsv_path, sep='\t', dtype={'Chr': str, 'Pos': int})
        gene_map['SiteID'] = gene_map['Chr'].astype(str) + ':' + gene_map['Pos'].astype(str)
        map_cols = ['SiteID', 'GeneSymbol', 'FeatureType']
        gene_map_unique = gene_map[map_cols].drop_duplicates(subset=['SiteID'])
        
        # Create a temporary SiteID column from the index for merging (Chr:Pos)
        temp_df = df.reset_index()
        temp_df['SiteID_key'] = temp_df['SiteID'].apply(lambda x: ':'.join(x.split(':')[0:2]))
        
        # Merge the gene/feature information
        temp_df = temp_df.merge(gene_map_unique, left_on='SiteID_key', right_on='SiteID', how='left', suffixes=(None, '_map'))
        
        temp_df.rename(columns={'FeatureType': 'Functional_Region', 'GeneSymbol': 'Gene'}, inplace=True)
        
        temp_df['Functional_Region'] = temp_df['Functional_Region'].fillna('Intergenic')
        temp_df['Gene'] = temp_df['Gene'].fillna('Intergenic')
        
        # Select final annotation columns and merge back to the original index
        final_annot_cols = ['SiteID', 'Functional_Region', 'Gene']
        gene_annot_df = temp_df[final_annot_cols].drop_duplicates(subset=['SiteID']).set_index('SiteID')
        
        # Merge the annotation back into the input dataframe (retaining quantification columns)
        result_df = pd.merge(df.reset_index(), gene_annot_df.reset_index(), on='SiteID', how='left').set_index('SiteID')

    except Exception as e:
        print(f"  Warning: Failed to load Gene Map: {e}", file=sys.stderr)
        # Fallback columns if merge fails
        result_df = df.copy()
        result_df['Functional_Region'] = 'Unknown'
        result_df['Gene'] = 'Unknown'

    return result_df

def rank_inverse_normal_transform(series):
    """Performs Rank-Inverse Normal Transformation (RINT) on a pandas Series."""
    if series.nunique() < 2 or series.std() == 0:
        return series
    
    rank = series.rank(method='average')
    n = len(series)
    norm_rank = (rank - 0.5) / n
    
    transformed = pd.Series(
        data=scipy.stats.norm.ppf(norm_rank), 
        index=series.index
    )
    
    return (transformed - transformed.mean()) / transformed.std()


# --- 3. Core Processing Logic ---

def load_and_aggregate_raw_calls(root_search_dir, individual_id, min_edit_level):
    """Loads, filters, aggregates, and finds consensus for one individual."""
    
    search_path = os.path.join(root_search_dir, individual_id, '**', f'{individual_id}_*_raw.tsv')
    all_files = glob.glob(search_path, recursive=True)
    
    if not all_files:
        raise FileNotFoundError(f"No raw call files found for individual {individual_id} in {root_search_dir}")

    raw_df_list = []
    
    for f in all_files:
        tool = 'REDML' if '_redml_raw.tsv' in f else 'REDItools'
        
        try:
            if tool == 'REDML':
                cols = ['Chr', 'Pos', 'Ref', 'Alt', 'VariantReads', 'TotalReads']
                df = pd.read_csv(f, sep='\t', usecols=cols, dtype={'Chr': str, 'Pos': int})
                df['EditLevel'] = df['VariantReads'] / df['TotalReads'].replace(0, pd.NA) 
                df = df.dropna(subset=['EditLevel']) 
            else: # REDItools
                cols = ['Region', 'Position', 'Reference', 'AllSubs', 'Frequency']
                df = pd.read_csv(f, sep='\t', usecols=cols, dtype={'Region': str, 'Position': int})
                df.rename(columns={'Region': 'Chr', 'Position': 'Pos'}, inplace=True)
                df = explode_reditools_substitutions(df.rename(columns={'Reference': 'Ref', 'Frequency': 'EditLevel'}))
                if df.empty: continue
            
            # Apply common filters
            df = df[df['EditLevel'] >= min_edit_level].copy()
            df = df[((df['Ref'] == 'A') & (df['Alt'] == 'G')) | ((df['Ref'] == 'T') & (df['Alt'] == 'C'))].copy()
            
            # Add metadata
            df['Tool'] = tool
            cell_type = os.path.basename(f).split(f'{individual_id}_')[1].split(f'_{tool.lower()}_raw.tsv')[0]
            df['CellType'] = cell_type
            df['SiteID'] = df['Chr'].astype(str) + ':' + df['Pos'].astype(str) + ':' + df['Ref'] + '>' + df['Alt']
            
            raw_df_list.append(df.loc[:, ['SiteID', 'CellType', 'Tool', 'EditLevel']])
            
        except Exception as e:
            print(f"Warning: Failed to process file {os.path.basename(f)}: {e}", file=sys.stderr)
            
    if not raw_df_list:
        raise ValueError(f"No valid data loaded after filtering for individual {individual_id}.")

    return pd.concat(raw_df_list, ignore_index=True)


def run_individual_processing(args):
    """Runs the full aggregation, annotation, and normalization pipeline for one individual."""

    raw_df = load_and_aggregate_raw_calls(args.root_search_dir, args.individual_id, args.min_edit_level)
    print(f"Total raw entries loaded for {args.individual_id} after Edit Level >= {args.min_edit_level} filter: {len(raw_df)}")

    # 1. Consensus Filter (RED-ML intersect REDItools) - PER-INDIVIDUAL
    consensus_check = raw_df.groupby('SiteID')['Tool'].nunique().reset_index(name='ToolCount')
    consensus_sites = consensus_check[consensus_check['ToolCount'] == 2]['SiteID']
    
    consensus_df = raw_df[raw_df['SiteID'].isin(consensus_sites)].copy()
    if consensus_df.empty:
        print("No consensus sites found for this individual. Skipping output.", file=sys.stderr)
        return
        
    print(f"Unique consensus sites for {args.individual_id}: {consensus_df['SiteID'].nunique()}")
    
    
    # 2. Aggregate and Pivot (Site x CellType Matrix)
    matrix_df = consensus_df.pivot(index='SiteID', columns='CellType', values='EditLevel')
    matrix_df.fillna(0, inplace=True)
    
    
    # 3. Annotation (CRITICAL: Using REDIPortal)
    print("Starting annotation...")
    annotated_df = matrix_df.copy()
    
    annotated_df = annotate_rediportal_status(annotated_df, args.rediportal_bed)
    
    # Functional Region and Gene Annotation
    final_annotated_df = annotate_functional_region(annotated_df, args.gene_map_tsv)

    annotation_cols = ['REDIPortal_Status', 'Functional_Region', 'Gene']
    quant_cols = [col for col in final_annotated_df.columns if col not in annotation_cols]
    quant_matrix = final_annotated_df[quant_cols]

    
    # 4. Calculate Gene Median and Normalize
    print("Calculating median gene editing ratio and normalizing...")
    
    # 4a. Calculate Median Editing Ratio per Gene
    median_group_df = final_annotated_df.reset_index()
    gene_median_df = median_group_df.groupby('Gene')[quant_cols].median()
    gene_median_df = gene_median_df[gene_median_df.index != 'Intergenic'].copy()
    
    # 4b. Normalize for FastQTL (RINT)
    normalized_quant_matrix = quant_matrix.apply(rank_inverse_normal_transform, axis=0)
    normalized_gene_median = gene_median_df.apply(rank_inverse_normal_transform, axis=0)

    
    # 5. Final Output
    print(f"Saving final matrix to {args.output_file}")
    with open(args.output_file, 'w') as f:
        f.write(f"# Processed Individual: {args.individual_id}\n")
        
        # Site-level matrix with annotation (Normalized)
        f.write("# --- SITE-LEVEL QUANTIFICATION (Normalized) ---\n")
        final_site_df = final_annotated_df[annotation_cols].merge(normalized_quant_matrix, 
                                                            left_index=True, right_index=True)
        final_site_df.to_csv(f, sep='\t', index=True)
        
        # Gene-level matrix (Normalized)
        f.write("\n# --- GENE-LEVEL MEDIAN QUANTIFICATION (Normalized for FastQTL) ---\n")
        normalized_gene_median.to_csv(f, sep='\t', index=True)

# --- 4. Main Execution ---
if __name__ == "__main__":
    try:
        run_individual_processing(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 3 Processing: {e}", file=sys.stderr)
        sys.exit(1)
