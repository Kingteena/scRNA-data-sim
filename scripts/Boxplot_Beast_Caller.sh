#!/bin/bash 
# This script calls upon Boxplot_beast.R with key input arguments; TRUETREE_PATH,BEASTTREE_PATH,OUTPUT_PATH and RUN_NAME
# For each run, the Origin CellCoal tree beign compared varies, but this is a feature of our experimentation due to time constraints of the project 
# Below lists which .vcf file and thus cellcoal tree file, each respective run compares its 100 simulated BEAST trees to; 

## RUN 1 - 0001.vcf / 0001.tree (2609 SNVs) 
## RUN 2 - 0003.vcf / 0003.tree (2608 SNVs) 

## RUN 3 - 0002.vcf / 0002.tree (2168 SNVs)
## RUN 4 - 0001.vcf / 0001.tree (2216 SNVs) 

## RUN 5 - 0001.vcf / 0001.tree (1726 SNVs) 
## RUN 6 - 0001.vcf / 0001.tree (1658 SNVs)

## RUN 7 - 0003.vcf / 0003.tree (1111 SNVs) 
## RUN 8 - 0003.vcf / 0003.tree (1028 SNVs) 


OUTPUT_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/beast
BOXPLOT_BEAST_SCRIPT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/boxplot_beast.R
 

for folder in "$OUTPUT_FOLDER"/run_3; do
        RUN_NUM=$(basename "$folder")
        #RUN_FOLDER=$OUTPUT_FOLDER/$RUN_NAME

        # Creating Output directory 
        mkdir -p "$folder/metric_calcs"

         #case statment for the exact cellcoal trees stated above  
        # Grouping specific runs as they share the same origin tree path
        case "$RUN_NUM" in
        "run_1" | "run_4" | "run_5" | "run_6")
            true_tree_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/cellcoal/${RUN_NUM}/trees_dir/0001.tree" 
            ;;
        "run_3")
            true_tree_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/cellcoal/${RUN_NUM}/trees_dir/0002.tree"
            ;;
        *)
            # Default path for all other runs (run_2, 7, 8)
            true_tree_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/cellcoal/${RUN_NUM}/trees_dir/0003.tree"
            ;;
    esac

    beast_trees_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/beast/${RUN_NUM}/tree"
    output_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/beast/${RUN_NUM}/metric_calcs"

    echo "Running analysis for $RUN_NUM..."

    Rscript $BOXPLOT_BEAST_SCRIPT "$true_tree_path" "$beast_trees_path" "$output_path" "$RUN_NUM"
done

