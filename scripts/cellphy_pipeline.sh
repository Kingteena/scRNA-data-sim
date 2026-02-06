#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 14
#SBATCH --mem=300
#SBATCH --time=02:00:00
#SBATCH --job-name=cellphy_array
#SBATCH --output=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellphy/log/%a.log
#SBATCH --array=1-100             # 100 VCF files per run × 8 runs

# --- CONFIGURATION ---
# Define your directories explicitly
BASE_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation"
CELLCOAL_BASE="${BASE_DIR}/output/cellcoal"
OUTPUT_BASE="${BASE_DIR}/output/cellphy"
CELLPHY_BIN="${BASE_DIR}/tools/cellphy-0.9.4/cellphy.sh"

# --- ARRAY MATH ---
# With 100 VCF files per run and 8 runs = 800 total tasks
# Calculate run number and VCF number directly from array task ID
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID / 100 ))
VCF_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

# Construct the VCF file path
VCF_FILE="${CELLCOAL_BASE}/run_${RUN_NUMBER}/expressed_snvs/${VCF_NAME}.vcf"

# Check if file exists (sanity check)
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found: $VCF_FILE"
    exit 1
fi

# --- SETUP DIRECTORIES ---
TARGET_OUTPUT_DIR="${OUTPUT_BASE}/run_${RUN_NUMBER}/${VCF_NAME}"
mkdir -p "$TARGET_OUTPUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: Processing run_${RUN_NUMBER} - ${VCF_NAME}"
echo "Input: $VCF_FILE"
echo "Output: $TARGET_OUTPUT_DIR"

# --- RUN CELLPHY ---
"$CELLPHY_BIN" -t 14 -o "$TARGET_OUTPUT_DIR" "$VCF_FILE"

echo "Finished with CellPhy run for run_${RUN_NUMBER}/${VCF_NAME}"
