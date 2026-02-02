#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=01:30:00
#SBATCH --job-name=beast_array
#SBATCH --output=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast/logs/%a.log  # %A is JobID, %a is ArrayIdx
#SBATCH --array=0-11              # 6 VCFs * 2 Runs = 12 Tasks

# --- CONFIGURATION ---
# Define your directories explicitly
BASE_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation"
SCRIPT_DIR="${BASE_DIR}/scripts/templates"
VCF_DIR="${BASE_DIR}/output/cellcoal/expressed_snvs"            
OUTPUT_BASE="${BASE_DIR}/output/beast"
PYTHON_SCRIPT="${SCRIPT_DIR}/vcf2xml.py"
XML_TEMPLATE="${SCRIPT_DIR}/test_gt16_error.xml"

NUM_REPLICATES=2  # Number of BEAST runs per VCF

# Load Modules
module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
module load beagle-lib/4.0.0-GCC-11.2.0

# --- ARRAY MATH ---
# 1. Get the list of VCF files
FILES=(${VCF_DIR}/*.vcf)

# 2. Determine which VCF to use (Index / 4)
# 0-3 -> File 0; 4-7 -> File 1; etc.
FILE_INDEX=$((SLURM_ARRAY_TASK_ID / $NUM_REPLICATES))
VCF_FILE=${FILES[$FILE_INDEX]}
FILE_NAME=$(basename "${VCF_FILE%.vcf}")

# 3. Determine which Run Number this is (Index % 4 + 1)
# 0 -> 1; 1 -> 2; 2 -> 3; 3 -> 4; 4 -> 1...
RUN_NUM=$((SLURM_ARRAY_TASK_ID % $NUM_REPLICATES + 1))

# Check if file exists (sanity check)
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found at index $FILE_INDEX"
    exit 1
fi

# --- SETUP DIRECTORIES ---
TARGET_XML_DIR="${OUTPUT_BASE}/${FILE_NAME}/xml"
TARGET_RUN_DIR="${OUTPUT_BASE}/${FILE_NAME}/run"
TARGET_TREE_DIR="${OUTPUT_BASE}/${FILE_NAME}/tree" # If your XML uses this, ensure it exists

mkdir -p "$TARGET_XML_DIR" "$TARGET_RUN_DIR" "$TARGET_TREE_DIR"

echo "Task $SLURM_ARRAY_TASK_ID: Processing $FILE_NAME | Run #$RUN_NUM"

# --- GENERATE XML ---
# We generate a unique temp XML for THIS task to avoid race conditions
TEMP_XML="${TARGET_XML_DIR}/temp_${RUN_NUM}.xml"
FINAL_XML="${TARGET_XML_DIR}/${RUN_NUM}.xml"

# 1. Run Python to create the base XML (Writing to a temp file unique to this run)
python3 "$PYTHON_SCRIPT" \
    "$VCF_FILE" \
    "$XML_TEMPLATE" \
    "$TEMP_XML" \
    --filtering_density 0.3

# 2. Run SED to inject the seed/path path
# Replaces $(seed) with the specific run directory for this instance
sed "s|\$(seed)|${TARGET_RUN_DIR}/${RUN_NUM}|g" "$TEMP_XML" > "$FINAL_XML"

# 3. Cleanup temp file
rm "$TEMP_XML"

# --- RUN BEAST ---
beast -beagle_CPU -beagle_SSE -overwrite "$FINAL_XML"