#!/bin/bash 
#need to add sbatch command here before we commit it to git 
#sbatch -N 1 -c 1 -n 1 --time=00:10:00 --mem=3GB -o logs/Cell_States.log
 
OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal
SCRIPT_RUN_EXPRESSIONS=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/run_expressions.sh

folder=$OUTPUTS_FOLDER
# for folder in $OUTPUTS_FOLDER/*; do
    echo "Processing folder: $(basename $folder)" 
    mkdir "$folder/expressed_snvs"
    for tree_file in $folder/tree_files/* ; do
        file_base=${tree_file##*.}
        vcf_file="$folder/cleaned_reps/${file_base}.vcf"
        vcf_output_file="$folder/expressed_snvs/expressed_snvs_${file_base}.vcf"
        echo "Calling expression simulation for tree: $tree_file and vcf: $vcf_file"
        bash $SCRIPT_RUN_EXPRESSIONS "$tree_file" "$vcf_file" "$vcf_output_file"
    done
# done

#folder=$OUTPUTS_FOLDER
# for folder in $OUTPUTS_FOLDER/*; do
    #echo "Processing folder: $(basename $folder)" 
    #mkdir "$folder/expressed_snvs"
    #for tree_file in $folder/tree_files/* ; do
        #file_base=$(basename "$tree_file" | cut -d. -f2)
        #vcf_file="$folder/cleaned_reps/${file_base}.vcf"
        #vcf_output_file="$folder/expressed_snvs/expressed_snvs_${file_base}.vcf"
        #echo "Calling expression simulation for tree: $tree_file and vcf: $vcf_file"
        #bash $SCRIPT_RUN_EXPRESSIONS "$tree_file" "$vcf_file" "$vcf_output_file"
    #done
# done


