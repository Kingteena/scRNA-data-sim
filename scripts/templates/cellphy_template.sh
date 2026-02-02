#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --time=01:30:00
#SBATCH --mem=2GB
#SBATCH -o logs/EstOutput.log

# No notifications given the fact that we don't want to spam users 

#Running CellPhy with the merged.vcf.gz file to find "best tree" with Maximum Likelihood Algorithm 

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf_file> <output_directory>"
    exit 1
fi

/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellphy-0.9.4/cellphy.sh -t 24 -o $2/CellPhy_Output $1