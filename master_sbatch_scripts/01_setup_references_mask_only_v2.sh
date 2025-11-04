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

# --- Directory Setup ---
BASE_DIR="/path/to/your/reference_data"
ANNOTATION_DIR="${BASE_DIR}/Annotation"
mkdir -p ${ANNOTATION_DIR}
mkdir -p logs

echo "Focusing on downloading and preparing repetitive element mask..."

# --- 1. Verify Existing Core Files (Optional Sanity Check) ---
# NOTE: This FASTA file MUST be the Ensembl-style (e.g., "1" not "chr1") for BAM compatibility.
GENOME_FASTA="${BASE_DIR}/Genome/GRCh38.no_alt.fa"
if [ ! -f "${GENOME_FASTA}" ]; then
    echo "WARNING: FASTA file not found at ${GENOME_FASTA}. Required for later stages."
    # If the FASTA is missing, you must uncomment and run the download/indexing steps from the original Phase 1.
fi

# --- 2. Download and Prepare Repetitive Element Mask (Blacklist) ---

echo -e "\nDownloading RepeatMasker elements from UCSC (GRCh38)..."
UCSC_RMSK_FILE="${ANNOTATION_DIR}/hg38_rmsk.txt.gz"
# Using wget -qO to quietly output to file. set -e ensures failure halts script.
wget -qO ${UCSC_RMSK_FILE} "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

# The check for download success is now implicitly handled by 'set -e' if wget fails.
# However, explicit check for CRITICAL failure remains good practice:
if [ ! -s "${UCSC_RMSK_FILE}" ]; then
    echo "CRITICAL ERROR: Downloaded file ${UCSC_RMSK_FILE} is missing or empty. Exiting."
    exit 1
fi

# Create the specific blacklist BED file (Simple Repeats and Alu elements)
BLACKLIST_BED="${ANNOTATION_DIR}/GRCh38_SimpleAlu_Blacklist.bed"

echo "Creating Simple Repeats and Alu Blacklist and converting to Ensembl chromosome names..."
# Pipeline: gunzip -> awk (filter/extract) -> sed (convert chr) -> sort (cleanup/standardize) -> file
gunzip -c ${UCSC_RMSK_FILE} | \
awk 'BEGIN {OFS="\t"} ($12 == "Simple_repeat" || $12 ~ /^Alu/) {print $6, $7-1, $8, $10}' | \
# --- THE CRITICAL FIX: Strip "chr" prefix for Ensembl BAM compatibility ---
sed 's/^chr//' | \
# --------------------------------------------------------------------------
sort -k1,1 -k2,2n > ${BLACKLIST_BED}

echo "Blacklist created with $(wc -l ${BLACKLIST_BED} | awk '{print $1}') regions."

echo -e "\nRevised Phase 1 Setup complete. Blacklist available at: ${BLACKLIST_BED}"
