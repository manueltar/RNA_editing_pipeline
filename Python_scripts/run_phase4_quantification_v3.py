#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob
import sys
import pysam
from collections import defaultdict
import re
import gzip
import numpy as np

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 4 (Reworked) v3: Per-Cell Quantification using Phase 3 Matrix with optional annotation filtering.")
parser.add_argument("--phase3_matrix", required=True, help="Input matrix from Phase 3 (INDIVIDUAL_ID_annotated_raw_matrix.tsv).")
parser.add_argument("--output_file", required=True, help="Path to save the final cell-type quantification matrix.")
parser.add_argument("--input_bams_dir", required=True, help="Directory containing the 31 input BAM files.")
parser.add_argument("--bam_pattern", default="scRNA_file_*.bam", help="Pattern to identify the 31 BAM files.")
parser.add_argument("--germline_vcf", required=True, help="Path to the individual's germline VCF file.")
parser.add_argument("--individual_id", required=True, help="Individual ID for VCF column query (e.g., IndA).")
parser.add_argument("--gtf_annotation", required=True, help="Path to the Ensembl GTF file.")
parser.add_argument("--splice_site_threshold", type=int, default=4, help="Exclude sites within this many bp of a splice junction.")
parser.add_argument("--min_read_coverage", type=int, default=10, help="Minimum TotalReads required at a site.")
parser.add_argument("--threads", type=int, default=16, help="Number of threads for parallel processing (optional).")
# NEW OPTIONAL FILTER ARGUMENTS
parser.add_argument("--filter_redip_status", type=str, default=None, help="Filter Phase 3 sites by a specific REDIPortal_Status (e.g., 'Known'). Set to None to disable.")
parser.add_argument("--filter_func_region", type=str, default=None, help="Filter Phase 3 sites by a specific Functional_Region (e.g., 'UTR3'). Set to None to disable.")
args = parser.parse_args()

BAM_FILES = sorted(glob.glob(os.path.join(args.input_bams_dir, args.bam_pattern)))
if not BAM_FILES:
    print("FATAL ERROR: No BAM files found based on pattern.", file=sys.stderr)
    sys.exit(1)

# --- 2. Filtering and Annotation Logic (Pysam/GTF Helpers) ---

# Global variable for GTF junctions
GTF_JUNCTION_MAP = None 

# --- Helper functions (parse_gtf_attributes, load_gtf_splice_junctions, 
# --- parse_gtf_and_get_annotation, check_germline_status, quantify_site_per_bam) 
# --- remain structurally the same as v2 but are included here for completeness.

def parse_gtf_attributes(attribute_str: str) -> dict:
    """Parses the attribute column of a GTF file."""
    attributes = {}
    matches = re.findall(r'\s*([^;]+?)\s+\"([^;]+?)\"\s*', attribute_str)
    for key, value in matches:
        attributes[key.strip()] = value.strip()
    return attributes

def load_gtf_splice_junctions(gtf_path: str) -> dict:
    """Loads all known splice junction coordinates (exon ends/starts) from GTF."""
    print("Loading GTF splice junction coordinates...")
    junctions = defaultdict(list)
    opener = gzip.open if gtf_path.endswith('.gz') else open
    
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] == 'exon':
                chrom = parts[0].replace('chr', '')
                start = int(parts[3])
                end = int(parts[4])
                junctions[chrom].append(start)
                junctions[chrom].append(end) 

    for chrom in junctions:
        junctions[chrom] = sorted(list(set(junctions[chrom])))
    return junctions

def parse_gtf_and_get_annotation(gtf_file, chrom, pos, splice_threshold):
    """Checks for proximity to splice sites and returns global filter status and annotation."""
    global GTF_JUNCTION_MAP
    if GTF_JUNCTION_MAP is None:
        try:
            GTF_JUNCTION_MAP = load_gtf_splice_junctions(gtf_file)
        except Exception as e:
            print(f"CRITICAL: Failed to load GTF for SJ filtering: {e}", file=sys.stderr)
            return {'GlobalFilterStatus': 'GTF_Error', 'MinDistToSplice': 9999} 

    chrom = chrom.replace('chr', '')
    min_dist = 9999
    if chrom in GTF_JUNCTION_MAP:
        junctions = np.array(GTF_JUNCTION_MAP[chrom])
        min_dist = np.min(np.abs(junctions - pos))
        
        if min_dist <= splice_threshold:
            return {
                'GlobalFilterStatus': f"SJ_Filtered_<{splice_threshold}bp", 
                'MinDistToSplice': int(min_dist)
            }
            
    return {
        'GlobalFilterStatus': 'PASS', 
        'MinDistToSplice': int(min_dist), 
        'FunctionalRegion': 'N/A', 
        'gene_name': 'N/A', 
        'transcript_id': 'N/A', 
        'gene_biotype': 'N/A'
    }

def check_germline_status(chrom, pos, ref, alt, vcf_path, ind_id):
    """Uses Pysam to query VCF for germline SNP status."""
    is_germline_snp = False
    vcf_status = "SomaticEdit"
    
    try:
        vcf_in = pysam.VariantFile(vcf_path)
        for rec in vcf_in.fetch(chrom, pos - 1, pos):
            if rec.ref == ref and alt in rec.alts:
                sample = rec.samples[ind_id]
                if sample["GT"] in [(0, 1), (1, 1), (1, 0)]: 
                    is_germline_snp = True
                    vcf_status = "GermlineSNP"
                    break
        vcf_in.close()

    except Exception as e:
        vcf_status = f"VCF_Error: {e}"
        print(f"Warning: Failed VCF query at {chrom}:{pos}: {e}", file=sys.stderr)

    return is_germline_snp, vcf_status

def quantify_site_per_bam(chrom, pos, ref, alt, bam_path, min_coverage):
    """Pysam pileup and count bases, applying the MIN_COVERAGE mask."""
    total_reads = 0
    variant_reads = 0
    
    try:
        bam_file = pysam.AlignmentFile(bam_path, "rb")
        for pileupcolumn in bam_file.pileup(chrom, pos - 1, pos, truncate=True, max_depth=100000):
            if pileupcolumn.pos == pos - 1:
                total_reads = pileupcolumn.n
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.query_position is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == alt:
                            variant_reads += 1
                break
        bam_file.close()

    except Exception as e:
        print(f"Error querying BAM {os.path.basename(bam_path)} at {chrom}:{pos}: {e}", file=sys.stderr)
        return {'TotalReads': 0, 'VariantReads': 0, 'EditingRatio': 'NA'}, "BAM_Error"

    if total_reads < min_coverage:
        return {
            'TotalReads': total_reads, 
            'VariantReads': variant_reads,
            'EditingRatio': 'NA' 
        }, "LowCoverage"
        
    editing_ratio = variant_reads / total_reads
    
    return {
        'TotalReads': total_reads, 
        'VariantReads': variant_reads,
        'EditingRatio': editing_ratio
    }, "PASS"


# --- 3. Main Execution and Matrix Construction ---

def run_phase4_quantification(args):
    """Loads Phase 3 matrix, applies **optional annotation filters**, applies global QC, then quantifies per-cell type."""
    
    print(f"Loading consensus sites from Phase 3 matrix: {args.phase3_matrix}")
    phase3_df = pd.read_csv(args.phase3_matrix, sep='\t', index_col='SiteID', comment='#')
    
    # --- STEP 1: Apply Optional Annotation Filters ---
    initial_site_count = len(phase3_df)
    
    if args.filter_redip_status and args.filter_redip_status.lower() != 'none':
        print(f"Filtering by REDIPortal_Status == '{args.filter_redip_status}'")
        phase3_df = phase3_df[phase3_df['REDIPortal_Status'] == args.filter_redip_status].copy()
        
    if args.filter_func_region and args.filter_func_region.lower() != 'none':
        print(f"Filtering by Functional_Region == '{args.filter_func_region}'")
        phase3_df = phase3_df[phase3_df['Functional_Region'] == args.filter_func_region].copy()
        
    final_site_count = len(phase3_df)
    print(f"Annotation Filter Summary: Sites reduced from {initial_site_count} to {final_site_count}.")
    
    if final_site_count == 0:
        print("No sites remaining after annotation filtering. Exiting gracefully.")
        pd.DataFrame().to_csv(args.output_file, sep='\t', index=False)
        return


    # --- STEP 2: Extract Coordinates and Apply Global QC (VCF/SJ) ---
    
    # Extract site components from the index (e.g., 'chr1:10000:A>G')
    site_components = phase3_df.index.to_series().str.split(':', expand=True)
    site_components.columns = ['Chr', 'Pos_str', 'Ref_Alt']
    site_components[['Ref', 'Alt']] = site_components['Ref_Alt'].str.split('>', expand=True)
    site_components['Pos'] = site_components['Pos_str'].astype(int)
    phase3_df[['Chr', 'Pos', 'Ref', 'Alt']] = site_components[['Chr', 'Pos', 'Ref', 'Alt']]

    print("Applying global VCF and Splice Junction filters...")
    site_annotation = {}
    
    for site_id, site in phase3_df.iterrows():
        chrom, pos, ref, alt = site.Chr, site.Pos, site.Ref, site.Alt
        
        # 1. Germline SNP Filter
        is_snp, vcf_status = check_germline_status(chrom, pos, ref, alt, args.germline_vcf, args.individual_id)
        
        # 2. Splice Junction Filter
        gtf_annot = parse_gtf_and_get_annotation(args.gtf_annotation, chrom, pos, args.splice_site_threshold)
        
        # Determine final Global Status
        if is_snp:
            global_status = vcf_status 
        elif gtf_annot['GlobalFilterStatus'].startswith('SJ_Filtered'):
            global_status = gtf_annot['GlobalFilterStatus']
        else:
            global_status = "PASS"
            
        site_annotation[site_id] = {
            'GlobalFilterStatus': global_status,
            'VCF_Status': vcf_status,
            'MinDistToSplice': gtf_annot.get('MinDistToSplice', 9999),
            # Keep original annotation carried from Phase 3
            'Phase3_FunctionalRegion': site.Functional_Region, 
            'Phase3_Gene': site.Gene,
            'Phase3_REDIPortal_Status': site.REDIPortal_Status
        }
        
    # --- STEP 3: Per-Cell Quantification and Matrix Building ---
    
    standard_cols = ['SiteID', 'Chr', 'Pos', 'Ref', 'Alt', 'GlobalFilterStatus', 
                     'VCF_Status', 'MinDistToSplice', 'Phase3_FunctionalRegion', 'Phase3_Gene', 'Phase3_REDIPortal_Status']
    
    matrix_rows = [] 
    
    print("Starting per-cell-type quantification on globally PASS sites...")

    for site_id, annot in site_annotation.items():
        # Only process sites that passed global VCF and SJ filters
        if annot['GlobalFilterStatus'] != 'PASS':
            continue

        site_row = phase3_df.loc[site_id]
        chrom, pos, ref, alt = site_row.Chr, site_row.Pos, site_row.Ref, site_row.Alt

        row = {'SiteID': site_id, 'Chr': chrom, 'Pos': pos, 'Ref': ref, 'Alt': alt}
        row.update(annot)
        
        # Quantification
        for bam_path in BAM_FILES:
            bam_id = os.path.basename(bam_path).replace('.bam', '')
            cell_type_id = bam_id.split('_')[-1] 
            
            metrics, status = quantify_site_per_bam(chrom, pos, ref, alt, bam_path, args.min_read_coverage)
            
            row[f'{cell_type_id}_ER'] = metrics['EditingRatio'] 
            row[f'{cell_type_id}_TR'] = metrics['TotalReads'] 
            row[f'{cell_type_id}_QC'] = status                 

        matrix_rows.append(row)

    # Convert the list of rows to the final DataFrame
    final_df = pd.DataFrame(matrix_rows)
    
    # Final cleanup and output
    final_cols = standard_cols + sorted([col for col in final_df.columns if col not in standard_cols])
    
    final_df = final_df[[col for col in final_cols if col in final_df.columns]]
    final_df.set_index('SiteID', inplace=True) # Set SiteID as the index

    final_df.to_csv(args.output_file, sep='\t', index=True, na_rep='NA')
    print(f"\nPhase 4 quantification complete. Final matrix size: {len(final_df)} rows.")


# --- 4. Main Execution ---
if __name__ == "__main__":
    try:
        run_phase4_quantification(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 4 Quantification: {e}", file=sys.stderr)
        sys.exit(1)

