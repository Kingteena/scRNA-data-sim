"""
This script is designed to run the expression processing logic for a specified simulation run. 
It reads .vcf files from the input directory, processes them using the `process_single_file` function from the `expression` module, and saves the results to an output directory. The script is intended to be run as part of a batch job on an HPC cluster, allowing for efficient processing of multiple files in parallel.
"""
import os
import glob
import time
import argparse
from tqdm import tqdm

from expression import process_single_file
def parse_args():
    """
    Parse command-line arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description="Run expression processing for a specified simulation run."
    )
    parser.add_argument(
        "--run",
        required=True,
        help="Run identifier (e.g., 4 for run_4).",
    )
    parser.add_argument(
        "--base-root",
        default="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal",
        help="Base output root containing run_* directories.",
    )
    return parser.parse_args()


def main():
    """
    Main function to execute the expression processing workflow.

    Configuration:
    - base_dir: Constructed from base_root and run identifier.
    - input_dir: Directory containing input .vcf files.
    - output_dir: Directory to save processed .vcf files.
    Simulation Config:
    - GENOME_LENGTH: Length of the simulated genome.
    - NUM_STATES: Number of cell states to simulate.
    - NUM_GENES: Total number of genes in the simulation.
    - ALPHA_GES: Parameter controlling gene expression simulation.
    - MINLIB: Minimum library size for simulated cells.
    - MAXLIB: Maximum library size for simulated cells.
    """
    args = parse_args()

    # --- CONFIGURATION ---
    base_dir = os.path.join(args.base_root, f"run_{args.run}")
    input_dir = base_dir + "/vcf_dir"
    output_dir = base_dir + "/expressed_snvs/"

    # Simulation Config
    sim_config = {
        'GENOME_LENGTH': 500000,
        'NUM_STATES': 5,
        'NUM_GENES': 30000,
        'ALPHA_GES': 0.5,
        'MINLIB': 1000,
        'MAXLIB': 5000
    }

    # Setup
    os.makedirs(output_dir, exist_ok=True)
    vcf_files = glob.glob(os.path.join(input_dir, "*.vcf"))

    if not vcf_files:
        print(f"No .vcf files found in {input_dir}")
        return

    print(f"Found {len(vcf_files)} files. Starting job...")
    start_time = time.time()

    for vcf_file in tqdm(vcf_files, desc="Processing VCFs"):

        # 1. Determine Output Path
        fname = os.path.basename(vcf_file)
        # fname = fname.replace("vcf.", "") + ".vcf"
        output_path = os.path.join(output_dir, fname)

        if os.path.exists(output_path):
            # print(f"Output already exists for {fname}, skipping.")
            continue
   
        # 2. Run the Library Logic
        result = process_single_file(vcf_file, output_path, sim_config)

        # Optional: Print errors immediately if they happen
        if result.lower().startswith("error"):
            print(f"\nError:  {result}")

    # Report
    duration = time.time() - start_time
    print(f"\nBatch Complete in {duration:.2f}s")

if __name__ == "__main__":
    main()
