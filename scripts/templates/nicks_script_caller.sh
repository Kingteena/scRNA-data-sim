#!/bin/bash 
#need to add sbatch command here before we commit it to git 
#sbatch -N 1 -c 1 -n 1 --time=00:10:00 --mem=3GB -o logs/Cell_States.log
 
OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal
SCRIPT_RUN_EXPRESSIONS=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/run_expressions.sh

folder=$OUTPUTS_FOLDER
mkdir "$folder/expressed_snvs"


for i in $folder/tree_files/*.trees; do
    tree_file="$folder/tree_files/$(printf "%04d" $i)".trees
    file_base=${tree_file%.trees}
    vcf_file="$folder/cleaned_reps/${file_base}.vcf"
    vcf_output_file="$folder/expressed_snvs/${file_base}.vcf"
    echo "Calling expression simulation for tree: $tree_file and vcf: $vcf_file"
    bash $SCRIPT_RUN_EXPRESSIONS "$tree_file" "$vcf_file" "$vcf_output_file"
done