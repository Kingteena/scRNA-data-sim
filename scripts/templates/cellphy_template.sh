#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --time=00:50:00
#SBATCH --mem=4GB
#SBATCH -o logs/EstOutput.log

# No notifications given the fact that we don't want to spam users

#Running CellPhy with the merged.vcf.gz file to find "best tree" with Maximum Likelihood Algorithm 

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf_file> <output_directory>"
    exit 1
fi

/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellphy-0.9.4/cellphy.sh -t 8  -o $2/CellPhy_Output $1