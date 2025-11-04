#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess

# --- 1. Configuration and Setup (Wrapper Interface) ---
parser = argparse.ArgumentParser(description="Wrapper for running REDItools3 on single-cell BAM chunk list.")
# Arguments are named to match the shell script's established pipeline flow:
parser.add_argument("--input_list", required=True, help="Path to the list file containing single-cell BAMs for this chunk.")
parser.add_argument("--fasta", required=True, help="Path to the reference FASTA file.")
parser.add_argument("--output", required=True, help="Path to save the raw REDItools output TSV.")
parser.add_argument("--threads", type=int, default=1, help="Number of threads for REDItools.")
parser.add_argument("--min_qual", type=int, default=20, help="Minimum base/mapping quality for calling.")
# Receives the explicit BED file argument from the shell
parser.add_argument("--bed_file", required=True, help="Path to the blacklist BED file (Alu/Simple Repeats).")
args = parser.parse_args()

# --- 2. Core Execution Function (Deployment Logic) ---
def execute_reditools_tool(input_list_file, fasta, output_file, threads, min_qual, bed_file):
    """
    Executes the external REDItools3 program using subprocess, passing arguments directly.
    REVISED: Filters are relaxed/removed to defer to Phase 4.
    """
    REDITOOLS_EXECUTABLE = "/path/to/reditools/REDItools3.py"

    print("--- REDItools3 Tool Execution ---")
    print(f"BAM List: {input_list_file}")
    print(f"Bed File: {bed_file}")
    print(f"Threads: {threads}, Min Base/Map Quality: {min_qual}")
    print(f"Output: {output_file}")

    # Build the full command list
    subprocess_command = [
        "python", REDITOOLS_EXECUTABLE,

        '--input-list', input_list_file,

        '-f', fasta,
        '-o', output_file,
        '-t', str(threads),

        '-q', str(min_qual), # Base Quality (Retained)
        '-u', str(min_qual), # Mapping Quality (Retained)
        
        # --- REVISED PARAMETERS FOR DEFERRED QC ---
        '-m', '0',  # Min number of reads supporting a variant -> RELAXED to 0
        '-c', '0',  # Min total coverage -> RELAXED to 0
        
        '--bed-file', bed_file, # Alu/Repeat Blacklist (Retained)
        
        # '-n', # CRITICAL: REMOVED. This filter for known variants is deferred to Phase 4.
        
        # NOTE: -l (Min edit level) defaults to 0.05 in REDItools3, which is sufficient sensitivity.
    ]

    try:
        print("\nExecuting command: " + " ".join(subprocess_command))
        result = subprocess.run(subprocess_command, check=True, text=True, capture_output=True)
        
        # ... (Log output and error handling) ...

    except FileNotFoundError:
        print(f"\nFATAL ERROR: REDItools3 script not found. Check the path: {REDITOOLS_EXECUTABLE}", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nFATAL ERROR: REDItools3 execution failed with return code {e.returncode}.", file=sys.stderr)
        sys.exit(1)

    # Final check: Was the output file successfully created?
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        return True

    print("\nFATAL ERROR: REDItools3 ran, but output file was not created or is empty.", file=sys.stderr)
    return False


# --- 3. Main Execution ---
if __name__ == "__main__":
    try:
        if execute_reditools_tool(args.input_list, args.fasta, args.output, args.threads, args.min_qual, args.bed_file):
            print(f"\nREDItools3 wrapper finished successfully. Output saved to {args.output}")
        else:
            sys.exit(1)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in REDItools3 wrapper: {e}", file=sys.stderr)
        sys.exit(1)
