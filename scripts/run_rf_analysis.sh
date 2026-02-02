#!/bin/bash
# Script to run RF distance analysis on all BEAST output trees

# Base directories
BEAST_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast"
CELLCOAL_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/tree_files"
OUTPUT_DIR="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/similarity"
RF_SCRIPT="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/RF_dist/V_RF_script.R"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each numbered directory in BEAST output
for beast_folder in "$BEAST_DIR"/[0-9]*; do
    if [ -d "$beast_folder" ]; then
        # Extract the replicate number (e.g., 0001 from /path/to/0001)
        num=$(basename "$beast_folder")
        
        # Define file paths
        true_tree="${CELLCOAL_DIR}/trees.${num}"
        inferred_tree="${beast_folder}/tree/combined.tree"
        output_png="${OUTPUT_DIR}/${num}.png"
        
        # Check if both input files exist
        if [ -f "$true_tree" ] && [ -f "$inferred_tree" ]; then
            echo "Processing replicate ${num}..."
            echo "  True tree: $true_tree"
            echo "  Inferred tree: $inferred_tree"
            echo "  Output: $output_png"
            
            # Run the R script
            Rscript "$RF_SCRIPT" "$true_tree" "$inferred_tree" "$output_png"
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Successfully processed ${num}"
            else
                echo "  ✗ Error processing ${num}"
            fi
            echo ""
        else
            echo "Skipping replicate ${num}: Missing input files"
            [ ! -f "$true_tree" ] && echo "  - True tree not found: $true_tree"
            [ ! -f "$inferred_tree" ] && echo "  - Inferred tree not found: $inferred_tree"
            echo ""
        fi
    fi
done

echo "RF analysis complete! Results saved to: $OUTPUT_DIR"
