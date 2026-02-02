#!/bin/bash 
 
OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output
METRICS_SCRIPT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/Metrics_Script.R
BEAST_OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast
CELLPHY_OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal


folder=$OUTPUTS_FOLDER
    echo "Running Comparative Metrics Scipt for each replicates Phylogenetic Trees" 
    mkdir "$OUTPUTS_FOLDER/similarity_tests"
    for true_tree in "$folder"/cellcoal_sahana/tree_files/* ; do
        # 1. Remove the directory path (everything up to the last /)
        filename="${true_tree##*/}"
        # 2. Remove the extension (everything after the last .)
        file_base="${filename%.*}"
        #file_base=$(basename "$true_tree" .tree)
        beast_tree="$folder/beast/${file_base}/tree/${file_base}.tree"
        cellphy_tree="$folder/cellphy/${file_base}/${file_base}.vcf.raxml.bestTree"
        mkdir $folder/similarity_tests/${file_base}
        output_path="$folder/similarity_tests/${file_base}"
        echo "Calling inference accuracy checker for True tree: $file_base and Inferred trees from CellPhy and BEAST "
        Rscript $METRICS_SCRIPT "$true_tree" "$beast_tree" "$cellphy_tree" "$output_path"  
    done
