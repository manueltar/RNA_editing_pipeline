#!/bin/bash
#SBATCH --job-name=Ref_Mask
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=debug

# --- CRITICAL ROBUSTNESS IMPROVEMENT ---
# Exit immediately if a command exits with a non-zero status.
set -e -o pipefail

# --- Environment Setup ---
module purge
module load wget/1.21             # Tool for downloading files
module load samtools/1.18         # For creating FASTA index (if needed)
module load bedtools/2.30         # Added: Needed to merge/sort blacklists

# --- Directory Setup ---
BASE_DIR="/path/to/your/reference_data"
ANNOTATION_DIR="${BASE_DIR}/Annotation"
mkdir -p ${ANNOTATION_DIR}
mkdir -p logs

echo "Focusing on downloading and preparing repetitive element mask..."

# --- 1. Verify Existing Core Files (Optional Sanity Check) ---
GENOME_FASTA="${BASE_DIR}/Genome/GRCh38.no_alt.fa"
if [ ! -f "${GENOME_FASTA}" ]; then
    echo "WARNING: FASTA file not found at ${GENOME_FASTA}. Required for later stages."
fi

# --- 2. Download and Prepare Repetitive Element Mask (Blacklist) ---

echo -e "\nDownloading RepeatMasker elements from UCSC (GRCh38)..."
UCSC_RMSK_FILE="${ANNOTATION_DIR}/hg38_rmsk.txt.gz"
wget -qO ${UCSC_RMSK_FILE} "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

if [ ! -s "${UCSC_RMSK_FILE}" ]; then
    echo "CRITICAL ERROR: Downloaded file ${UCSC_RMSK_FILE} is missing or empty. Exiting."
    exit 1
fi

# --- CRITICAL FIX: Separate Simple Repeats and Alu Blacklists ---

# 2a. Create Simple Repeats Blacklist (for RED-ML --simpleRepeat)
SIMPLE_REPEATS_BED="${ANNOTATION_DIR}/01_REPEATS_Blacklist.bed"
echo "Creating Simple Repeats Blacklist: ${SIMPLE_REPEATS_BED}"

gunzip -c ${UCSC_RMSK_FILE} | \
awk 'BEGIN {OFS="\t"} ($12 == "Simple_repeat") {print $6, $7-1, $8, $10}' | \
# Strip "chr" prefix for Ensembl BAM compatibility
sed 's/^chr//' | \
sort -k1,1 -k2,2n > ${SIMPLE_REPEATS_BED}

# 2b. Create Alu Blacklist (for RED-ML --alu)
ALU_BED="${ANNOTATION_DIR}/01_ALU_Blacklist.bed"
echo "Creating Alu Blacklist: ${ALU_BED}"

gunzip -c ${UCSC_RMSK_FILE} | \
awk 'BEGIN {OFS="\t"} ($12 ~ /^Alu/) {print $6, $7-1, $8, $10}' | \
# Using /^Alu/ to capture the entire family.
# Strip "chr" prefix for Ensembl BAM compatibility
sed 's/^chr//' | \
sort -k1,1 -k2,2n > ${ALU_BED}

# 2c. Create a Merged Master Blacklist (for REDItools -r)
# REDItools often takes a single file for exclusion regions.
MASTER_BLACKLIST_BED="${ANNOTATION_DIR}/GRCh38_SimpleAlu_Master_Blacklist.bed"
echo "Creating combined Master Blacklist for REDItools..."

# Concatenate the two files, sort, and merge overlapping regions to create a clean mask
cat ${SIMPLE_REPEATS_BED} ${ALU_BED} | \
sort -k1,1 -k2,2n | \
bedtools merge -i - > ${MASTER_BLACKLIST_BED}

echo "Revised Phase 1 Setup complete."
echo "Simple Repeats: $(wc -l ${SIMPLE_REPEATS_BED} | awk '{print $1}') regions."
echo "Alu Repeats: $(wc -l ${ALU_BED} | awk '{print $1}') regions."
echo "Master Blacklist: $(wc -l ${MASTER_BLACKLIST_BED} | awk '{print $1}') merged regions."
