"""
Main Execution Script.
Runs the expression_lib functions in parallel.
"""
import os
import glob
import time
from tqdm import tqdm

# Import your custom library
from expression import process_single_file

def main():
    # --- CONFIGURATION ---
    BASE_DIR = "/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/run_4"
    INPUT_DIR = BASE_DIR + "/vcf_dir"
    OUTPUT_DIR = BASE_DIR + "/expressed_snvs/"
    
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
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    vcf_files = glob.glob(os.path.join(INPUT_DIR, "*.vcf"))
    
    if not vcf_files:
        print(f"No .vcf files found in {INPUT_DIR}")
        return

    print(f"Found {len(vcf_files)} files. Starting job...")
    start_time = time.time()

    for vcf_file in tqdm(vcf_files, desc="Processing VCFs"):
        
        # 1. Determine Output Path
        fname = os.path.basename(vcf_file)
        # fname = fname.replace("vcf.", "") + ".vcf"
        output_path = os.path.join(OUTPUT_DIR, fname)

        if os.path.exists(output_path):
            # print(f"Output already exists for {fname}, skipping.")
            continue
        
        # 2. Run the Library Logic
        result = process_single_file(vcf_file, output_path, sim_config)
        
        # Optional: Print errors immediately if they happen
        if "ERROR" in result:
            print(f"\nError:  {result}")

    # Report
    duration = time.time() - start_time
    print(f"\nBatch Complete in {duration:.2f}s")

if __name__ == "__main__":
    main()