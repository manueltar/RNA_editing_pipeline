#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob
import sys
import gzip
import re

# --- 1. Argument Parsing and Setup ---
parser = argparse.ArgumentParser(description="Phase 3: Individual Aggregation, Consensus Filtering, and Annotation.")
parser.add_argument("--individual_id", required=True, 
                    help="The unique ID of the individual being processed (e.g., IND_0001).")
parser.add_argument("--root_search_dir", required=True, 
                    help="Top-level directory to recursively search for raw call files.")
parser.add_argument("--output_file", required=True, 
                    help="Path to save the final annotated matrix with raw (unnormalized) Edit Levels for this individual.")
parser.add_argument("--min_edit_level", type=float, default=0.1, 
                    help="Minimum Editing Level (e.g., 0.1 = 10%%) required for a site to be kept.")
parser.add_argument("--rediportal_bed", required=True,
                    help="Path to the REDIPortal Known Sites BED file.")
parser.add_argument("--ensembl_gtf", required=True,
                    help="Path to the Ensembl GTF file (e.g., Homo_sapiens.GRCh38.109.gtf.gz).")
args = parser.parse_args()


# --- 2. Annotation Helper Functions ---

def parse_gtf_attributes(attribute_str: str) -> dict:
    """Parses the attribute column of a GTF file."""
    attributes = {}
    matches = re.findall(r'\s*([^;]+?)\s+\"([^;]+?)\"\s*', attribute_str)
    
    for key, value in matches:
        attributes[key.strip()] = value.strip()
    return attributes

def load_gtf_features(gtf_path: str) -> pd.DataFrame:
    """Loads and preprocesses the GTF file to create a feature map."""
    print("  -> Loading and parsing Ensembl GTF file. This may take time...")
    
    TARGET_FEATURES = {'exon', 'UTR', 'CDS'} 
    FEATURE_TYPES = {'exon': 'Exon', 'CDS': 'CDS', 'five_prime_utr': 'UTR5', 'three_prime_utr': 'UTR3'}

    gtf_records = []
    opener = gzip.open if gtf_path.endswith('.gz') else open
    
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
                
            feature_type = parts[2]
            
            if feature_type in TARGET_FEATURES or feature_type.endswith('_utr'):
                attributes = parse_gtf_attributes(parts[8])
                gene_symbol = attributes.get('gene_name', attributes.get('gene_id', 'Unknown'))

                gtf_records.append({
                    'Chr': parts[0].replace('chr', ''), 
                    'Start': int(parts[3]), 
                    'End': int(parts[4]),    
                    'Feature': FEATURE_TYPES.get(feature_type.lower(), feature_type),
                    'Strand': parts[6],
                    'Gene': gene_symbol
                })

    gtf_df = pd.DataFrame(gtf_records)
    gtf_df['Chr'] = gtf_df['Chr'].astype(str).str.replace('^chr', '', regex=True)
    gtf_df.sort_values(by=['Chr', 'Start'], inplace=True)
    
    print(f"  -> Loaded {len(gtf_df)} GTF features.")
    return gtf_df

# Global variable for lazy loading
GTF_FEATURE_MAP = None
def get_gtf_feature_map(gtf_path: str) -> pd.DataFrame:
    global GTF_FEATURE_MAP
    if GTF_FEATURE_MAP is None:
        GTF_FEATURE_MAP = load_gtf_features(gtf_path)
    return GTF_FEATURE_MAP

def annotate_rediportal_status(df: pd.DataFrame, rediportal_bed_path: str) -> pd.DataFrame:
    """Adds a 'REDIPortal_Status' column ('Known'/'Novel')."""
    df['REDIPortal_Status'] = 'Novel'
    try:
        rediportal_sites = pd.read_csv(rediportal_bed_path, sep='\t', header=None, usecols=[0,1], compression='infer', dtype={0: str})
        rediportal_sites.columns = ['Chr', 'Start']
        rediportal_sites['Pos'] = rediportal_sites['Start'] + 1 
        rediportal_sites['KnownID'] = rediportal_sites['Chr'].astype(str) + ':' + rediportal_sites['Pos'].astype(str)
        known_site_ids = set(rediportal_sites['KnownID'])
        
        df.loc[df.index.to_series().apply(lambda x: x.split(':')[0] + ':' + x.split(':')[1]).isin(known_site_ids), 'REDIPortal_Status'] = 'Known'
    except Exception as e:
        print(f"  Warning: Failed to load REDIPortal: {e}", file=sys.stderr)
    return df

def annotate_functional_region(df: pd.DataFrame, gtf_path: str) -> pd.DataFrame:
    """Annotates Functional_Region and Gene using GTF feature overlap."""
    print("  -> Annotating Functional Region and Gene using GTF...")
    gtf_df = get_gtf_feature_map(gtf_path)
    
    site_positions = df.index.to_series().apply(lambda x: x.split(':')).apply(pd.Series).rename(columns={0: 'Chr', 1: 'Pos_str', 2: 'Ref_Alt'})
    site_positions['Pos'] = site_positions['Pos_str'].astype(int)
    site_positions.set_index(df.index, inplace=True)

    df['Functional_Region'] = 'Intergenic'
    df['Gene'] = 'Intergenic'

    for site_id, row in site_positions.iterrows():
        chrom = row['Chr']
        pos = row['Pos']
        
        chrom_features = gtf_df[gtf_df['Chr'] == chrom]
        
        # Check for overlap: Start <= Pos <= End
        overlapping_features = chrom_features[
            (chrom_features['Start'] <= pos) & 
            (chrom_features['End'] >= pos)
        ]
        
        if not overlapping_features.empty:
            # Priority: CDS > UTR5/UTR3 > Exon
            
            cds_features = overlapping_features[overlapping_features['Feature'] == 'CDS']
            if not cds_features.empty:
                df.loc[site_id, 'Functional_Region'] = 'CDS'
                df.loc[site_id, 'Gene'] = cds_features['Gene'].iloc[0]
                continue
            
            utr_features = overlapping_features[overlapping_features['Feature'].isin(['UTR5', 'UTR3'])]
            if not utr_features.empty:
                best_feature = 'UTR3' if 'UTR3' in utr_features['Feature'].values else 'UTR5'
                utr_row = utr_features[utr_features['Feature'] == best_feature].iloc[0]
                df.loc[site_id, 'Functional_Region'] = best_feature
                df.loc[site_id, 'Gene'] = utr_row['Gene']
                continue
            
            exon_features = overlapping_features[overlapping_features['Feature'] == 'Exon']
            if not exon_features.empty:
                df.loc[site_id, 'Functional_Region'] = 'Exon'
                df.loc[site_id, 'Gene'] = exon_features['Gene'].iloc[0]
                continue

            # Fallback for unexpected feature types
            df.loc[site_id, 'Functional_Region'] = overlapping_features['Feature'].iloc[0]
            df.loc[site_id, 'Gene'] = overlapping_features['Gene'].iloc[0]

    return df


# --- 3. Core Processing Logic ---

def explode_reditools_substitutions(df: pd.DataFrame) -> pd.DataFrame:
    """Explodes the REDItools 'AllSubs' column into multiple rows, standardizing columns."""
    # (Same function as before, included for completeness)
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
            
            df = df[df['EditLevel'] >= min_edit_level].copy()
            df = df[((df['Ref'] == 'A') & (df['Alt'] == 'G')) | ((df['Ref'] == 'T') & (df['Alt'] == 'C'))].copy()
            
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
    """Runs the full aggregation, consensus filtering, and annotation pipeline."""

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
    
    
    # 2. Aggregate and Pivot (Site x CellType Matrix with RAW Edit Levels)
    matrix_df = consensus_df.pivot(index='SiteID', columns='CellType', values='EditLevel')
    matrix_df.fillna(0, inplace=True)
    
    
    # 3. Annotation
    print("Starting annotation...")
    annotated_df = matrix_df.copy()
    
    # REDIPortal Annotation
    annotated_df = annotate_rediportal_status(annotated_df, args.rediportal_bed)
    
    # Functional Region and Gene Annotation (Uses GTF)
    final_annotated_df = annotate_functional_region(annotated_df, args.ensembl_gtf)

    
    # 4. Final Output (Raw Matrix + Annotation)
    annotation_cols = ['REDIPortal_Status', 'Functional_Region', 'Gene']
    quant_cols = [col for col in final_annotated_df.columns if col not in annotation_cols]
    
    final_output_df = final_annotated_df[annotation_cols + quant_cols].copy()

    print(f"Saving final annotated raw matrix to {args.output_file}")
    
    # Write to file
    with open(args.output_file, 'w') as f:
        f.write(f"# Processed Individual: {args.individual_id}\n")
        f.write("# --- SITE-LEVEL QUANTIFICATION (Raw Edit Levels) ---\n")
        final_output_df.to_csv(f, sep='\t', index=True)


# --- 4. Main Execution ---
if __name__ == "__main__":
    try:
        run_individual_processing(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 3 Processing: {e}", file=sys.stderr)
        sys.exit(1)
