#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob
import sys
import pysam
from collections import defaultdict

# --- 1. Configuration and Setup ---
parser = argparse.ArgumentParser(description="Phase 4 (Reworked): Per-Cell Quantification, Germline, SJ Filter, and Annotation.")
parser.add_argument("--master_sites", required=True, help="Input file from Phase 3 (consensus_master_sites_p3.tsv).")
parser.add_argument("--output_file", required=True, help="Path to save the final cell-type quantification matrix.")
parser.add_argument("--input_bams_dir", required=True, help="Directory containing the 31 input BAM files.")
parser.add_argument("--bam_pattern", default="scRNA_file_*.bam", help="Pattern to identify the 31 BAM files.")
parser.add_argument("--germline_vcf", required=True, help="Path to the individual's germline VCF file (VCF filtering).")
parser.add_argument("--individual_id", required=True, help="Individual ID for VCF column query (e.g., IndA).")
parser.add_argument("--gtf_annotation", required=True, help="Path to the Ensembl GTF file (Annotation & SJ filter).")
parser.add_argument("--splice_site_threshold", type=int, default=4, help="Exclude sites within this many bp of a splice junction.")
parser.add_argument("--min_read_coverage", type=int, default=10, help="Minimum TotalReads required at a site to report a valid Editing Ratio (masking QC).")
parser.add_argument("--threads", type=int, default=16, help="Number of threads for parallel processing (optional).")
args = parser.parse_args()

# Global list of BAM files identified
BAM_FILES = sorted(glob.glob(os.path.join(args.input_bams_dir, args.bam_pattern)))
if not BAM_FILES:
    print("FATAL ERROR: No BAM files found based on pattern.", file=sys.stderr)
    sys.exit(1)

# --- 2. Filtering and Annotation Logic (Production Placeholders for VCF/GTF) ---

# NOTE: The check_germline_status and parse_gtf_and_get_annotation functions 
# remain the same as they do not depend on MIN_READ_COVERAGE.

def check_germline_status(chrom, pos, ref, alt, vcf_path, ind_id):
    """(Germline VCF check - Placeholder for pysam.VariantFile logic)"""
    # SIMULATION ONLY (Replace with real VCF lookup):
    if (pos % 20) == 0:
        return True, "GermlineSNP"
    return False, "SomaticEdit"
    # ... error handling logic ...

def parse_gtf_and_get_annotation(gtf_file, chrom, pos, splice_threshold):
    """(SJ Filter and GTF Annotation - Placeholder for GTF parser logic)"""
    # SIMULATION ONLY (Replace with real GTF lookup):
    dist = (pos % 25)
    if dist <= splice_threshold:
        return {'GlobalFilterStatus': f"SJ_Filtered (<={splice_threshold}bp)", 'MinDistToSplice': dist}
    region = 'Intronic'
    if 5 <= dist <= 10: region = 'three_prime_utr'
    elif 11 <= dist <= 20: region = 'CDS'
    elif 21 <= dist <= 24: region = 'five_prime_utr'
    return {
        'GlobalFilterStatus': 'PASS', 'MinDistToSplice': dist, 'FunctionalRegion': region,
        'gene_name': f'GENE_{chrom}', 'transcript_id': f'ENST_{pos}', 'gene_biotype': 'protein_coding'
    }

def quantify_site_per_bam(chrom, pos, ref, alt, bam_path, min_coverage):
    """
    PRODUCTION DEPLOYABLE LOGIC: Uses Pysam to perform pileup and count bases.
    Accepts min_coverage as a parameter.
    """
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


    # 1. Apply TotalReads >= MIN_COVERAGE Filter (Per-Cell QC)
    if total_reads < min_coverage:
        return {
            'TotalReads': total_reads, 
            'VariantReads': variant_reads,
            'EditingRatio': 'NA' # Mask the ratio
        }, "LowCoverage"
        
    # If passed coverage QC
    editing_ratio = variant_reads / total_reads
    
    return {
        'TotalReads': total_reads, 
        'VariantReads': variant_reads,
        'EditingRatio': editing_ratio
    }, "PASS"


# --- 3. Main Execution and Matrix Construction ---

def run_phase4_quantification(args):
    """Builds the final quantification matrix and applies all filters."""
    
    master_sites = pd.read_csv(args.master_sites, sep='\t')
    
    site_annotation = {}
    
    # --- Step 1: Global Filter Check (VCF & SJ) and Pre-calculate Annotation ---
    # ... (VCF/GTF global filter loop remains the same) ...

    # --- Step 2: Per-Cell Quantification and Matrix Building ---
    
    # Define standard site and annotation columns
    standard_cols = ['Chr', 'Pos', 'Ref', 'Alt', 'GlobalFilterStatus', 
                     'MinDistToSplice', 'FunctionalRegion', 
                     'gene_name', 'transcript_id', 'gene_biotype']
    
    matrix_rows = [] 

    for _, site in master_sites.iterrows():
        site_key = f"{site.Chr}:{site.Pos}:{site.Ref}>{site.Alt}"
        
        # Only process sites that passed global VCF and SJ filters
        if site_key not in site_annotation or site_annotation[site_key].get('GlobalFilterStatus') != 'PASS':
            continue

        row = site.to_dict()
        row.update(site_annotation[site_key])
        
        # Quantification for sites that passed global filters
        for bam_path in BAM_FILES:
            bam_id = os.path.basename(bam_path).replace('.bam', '')
            
            # KEY CHANGE: Pass the MIN_READ_COVERAGE parameter
            metrics, status = quantify_site_per_bam(site.Chr, site.Pos, site.Ref, site.Alt, bam_path, args.min_read_coverage)
            
            # Add metrics to the row
            row[f'{bam_id}_ER'] = metrics['EditingRatio']
            row[f'{bam_id}_TR'] = metrics['TotalReads']
        
        matrix_rows.append(row)

    # Convert the list of rows to the final DataFrame
    final_df = pd.DataFrame(matrix_rows)
    
    # Final cleanup and output
    final_cols = standard_cols + [col for col in final_df.columns if col not in standard_cols]
    
    final_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"\nPhase 4 quantification complete. Final matrix size: {len(final_df)} rows.")


# --- 4. Main Execution ---
if __name__ == "__main__":
    # Simplified execution block
    try:
        run_phase4_quantification(args)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in Phase 4 Quantification: {e}", file=sys.stderr)
        sys.exit(1)
