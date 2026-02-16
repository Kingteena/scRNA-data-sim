#!/bin/bash 
# This script calls upon boxplot_draw.R with key input arguments; true_trees_path, beast_trees_path, output_path and RUN_NUM
# For each run, the origin CellCoal tree for the respective 100 simulated datasets is compared to its corresponding Cellphy and BEAST inferred tree


OUTPUT_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output  #change to your own output folder 
BOXPLOT_BEAST_SCRIPT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/src/data_analysis/boxplot_draw.R  #change to directory of boxplot_draw.R script 
 

for folder in "$OUTPUT_FOLDER"/beast/run_7; do
        RUN_NUM=$(basename "$folder")
        RUN_FOLDER="$OUTPUT_FOLDER/$RUN_NUM"

        # Creating Output directory 
        OUTPUT_DIR="$OUTPUT_FOLDER/metrics_calcs/$RUN_NUM"
        echo "Creating output directory: $OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"

    true_trees_path="$OUTPUT_FOLDER/cellcoal/${RUN_NUM}/trees_dir"
    beast_trees_path="$OUTPUT_FOLDER/beast/${RUN_NUM}/trees"
    cellphy_trees_path="$OUTPUT_FOLDER/cellphy/${RUN_NUM}"
    output_path="$OUTPUT_FOLDER/metrics_calcs/${RUN_NUM}"

    echo "Running analysis for $RUN_NUM..."

    Rscript $BOXPLOT_BEAST_SCRIPT "$true_trees_path" "$beast_trees_path" "$cellphy_trees_path" "$output_path" "$RUN_NUM"
done

