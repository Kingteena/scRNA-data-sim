#!/bin/bash 
#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH -n 1 
#SBATCH --time=00:10:00 
#SBATCH --mem=3GB 
#SBATCH -o /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellphy/log/Cell_States.log
 
OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal
SCRIPT_RUN_EXPRESSIONS=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/run_expressions.sh
BEAST_OUTPUTS_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast

cd $(dirname $SCRIPT_RUN_EXPRESSIONS)

for folder in $OUTPUTS_FOLDER/*; do
    echo "Processing folder: $folder"
    mkdir -p "$folder/expressed_snvs"

    mkdir -p "$BEAST_OUTPUTS_FOLDER/$(basename $folder)/xml"
    mkdir -p "$BEAST_OUTPUTS_FOLDER/$(basename $folder)/run"
    mkdir -p "$BEAST_OUTPUTS_FOLDER/$(basename $folder)/tree"

    for tree_file in $folder/trees_dir/*.tree; do
        file_base=$(basename ${tree_file%.tree})
        vcf_file="$folder/cleaned_reps/${file_base}.vcf"
        vcf_output_file="$folder/expressed_snvs/${file_base}.vcf"
        echo "Calling expression simulation for tree: $tree_file and vcf: $vcf_file"
        bash $SCRIPT_RUN_EXPRESSIONS "$tree_file" "$vcf_file" "$vcf_output_file" "$BEAST_OUTPUTS_FOLDER/$(basename $folder)/xml/${file_base}.xml"
    done
done
