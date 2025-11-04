#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess

# --- 1. Configuration and Setup (Wrapper Interface) ---
parser = argparse.ArgumentParser(description="Wrapper for running RED-ML on single-cell BAM chunk list.")
# Arguments are named to match the shell script's established pipeline flow:
parser.add_argument("--input_list", required=True, help="Path to the list file containing single-cell BAMs for this chunk.")
parser.add_argument("--fasta", required=True, help="Path to the reference FASTA file.")
parser.add_argument("--output", required=True, help="Path to save the raw RED-ML output TSV.")
parser.add_argument("--threads", type=int, default=4, help="Number of threads for RED-ML.")
parser.add_argument("--min_qual", type=int, default=20, help="Minimum base quality for calling (Maps to RED-ML's --min_base_q).")
parser.add_argument("--blacklist", required=True, help="Path to the blacklist BED file (Alu/Simple Repeats).")
args = parser.parse_args()

# --- 2. Core Execution Function (Deployment Logic) ---
def execute_redml_tool(input_list_file, fasta, output_file, threads, min_qual, blacklist_bed):
    """
    Executes the external RED-ML tool using subprocess, mapping wrapper args to tool args.
    """
    # ASSUMED: The executable is the primary python script from the repo,
    # and it is accessible in the environment.
    REDML_EXECUTABLE = "python /path/to/redml/RED-ML.py" 

    print("--- RED-ML Tool Execution ---")
    print(f"BAM List: {input_list_file}")
    print(f"Blacklist: {blacklist_bed}")
    print(f"Threads: {threads}, Min Base Quality: {min_qual}")
    print(f"Output: {output_file}")

    # Build the full command list, mapping our wrapper arguments to the RED-ML tool's actual flags.
    # NOTE: --input-list and --blacklist are ASSUMED to be functional in your local RED-ML version.
    # The actual RED-ML tool arguments are used: -f, -t, --min_base_q.
    subprocess_command = [
        # Must pass the command and all arguments as separate strings
        "python", "/path/to/redml/RED-ML.py", 
        
        '--input-list', input_list_file,  # ASSUMED: Argument for chunk processing
        
        '-f', fasta, # Actual RED-ML argument
        '-o', output_file,
        '-t', str(threads), # Actual RED-ML argument
        
        '--min_base_q', str(min_qual), # Corrected RED-ML argument
        
        '--blacklist', blacklist_bed, # ASSUMED: Argument for Alu/Repeat filtering
        
        '--min_map_q', '20' # Recommended minimum mapping quality (common filter)
    ]

    try:
        print("\nExecuting command: " + " ".join(subprocess_command))
        # Execute the external program. 'check=True' ensures failure halts the script.
        result = subprocess.run(subprocess_command, check=True, text=True, capture_output=True)

        # Log output from the executable
        if result.stdout:
             print("\nRED-ML STDOUT:\n", result.stdout)
        if result.stderr:
             print("\nRED-ML STDERR:\n", result.stderr)

    except FileNotFoundError:
        print(f"\nFATAL ERROR: Python or RED-ML script not found. Check the path.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nFATAL ERROR: RED-ML execution failed with return code {e.returncode}.", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    # Final check: Was the output file successfully created?
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        return True
    
    print("\nFATAL ERROR: RED-ML ran, but output file was not created or is empty.", file=sys.stderr)
    return False


# --- 3. Main Execution ---
if __name__ == "__main__":
    try:
        if execute_redml_tool(args.input_list, args.fasta, args.output, args.threads, args.min_qual, args.blacklist):
            print(f"\nRED-ML wrapper finished successfully. Output saved to {args.output}")
        else:
            sys.exit(1)
    except Exception as e:
        print(f"\nFATAL UNCAUGHT ERROR in RED-ML wrapper: {e}", file=sys.stderr)
        sys.exit(1)
