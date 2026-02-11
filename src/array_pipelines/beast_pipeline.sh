#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=03:00:00
#SBATCH --job-name=beast_array
#SBATCH --output=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast/logs/%a.log
#SBATCH --array=73,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,300,301,302,303,304,305,306,307,317,318,329,330,331,332,333,334,344,345,346,347,348,355,356,357,358,359,360,361,369,370,371,372,375,376,377,378,382,383,384,385,386,387,388,389,399 
              # Adjust based on total number of VCF files



# --- CONFIGURATION ---
# Define your directories explicitly
BASE_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation"
SCRIPT_DIR="${BASE_DIR}/scripts/templates"
CELLCOAL_BASE="${BASE_DIR}/output/cellcoal"
OUTPUT_BASE="${BASE_DIR}/output/beast"
PYTHON_SCRIPT="${SCRIPT_DIR}/vcf2xml.py"
XML_TEMPLATE="${SCRIPT_DIR}/test_gt16_error.xml"

# Load Modules
module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
module load beagle-lib/4.0.0-GCC-11.2.0

# --- ARRAY MATH ---
# With 100 VCF files per run and 8 runs = 800 total tasks
# Calculate run number and VCF number directly from array task ID
RUN_NUMBER=$((SLURM_ARRAY_TASK_ID / 100 + 1))
VCF_INDEX=$((SLURM_ARRAY_TASK_ID % 100 + 1))
VCF_NAME=$(printf "%04d" $VCF_INDEX)

# Construct the VCF file path
VCF_FILE="${CELLCOAL_BASE}/run_${RUN_NUMBER}/expressed_snvs/${VCF_NAME}.vcf"

# Check if file exists (sanity check)
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found: $VCF_FILE"
    exit 1
fi


# --- SETUP DIRECTORIES ---
TARGET_XML_DIR="${OUTPUT_BASE}/run_${RUN_NUMBER}/xml"
TARGET_RUN_DIR="${OUTPUT_BASE}/run_${RUN_NUMBER}/run"

mkdir -p "$TARGET_XML_DIR" "$TARGET_RUN_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: Processing $RUN_NUMBER - $VCF_NAME"

# --- GENERATE XML ---
TEMP_XML="${TARGET_XML_DIR}/${VCF_NAME}_temp.xml"
FINAL_XML="${TARGET_XML_DIR}/${VCF_NAME}.xml"



# 1. Run Python to create the base XML
python3 "$PYTHON_SCRIPT" \
    "$VCF_FILE" \
    "$XML_TEMPLATE" \
    "$TEMP_XML" \
    --filtering_density 0.3

# 2. Run SED to inject the seed/path
sed "s|\$(seed)|${TARGET_RUN_DIR}/${VCF_NAME}|g" "$TEMP_XML" > "$FINAL_XML"

# 3. Cleanup temp file
rm "$TEMP_XML"

# --- RUN BEAST ---
beast -beagle_CPU -beagle_SSE -overwrite "$FINAL_XML"