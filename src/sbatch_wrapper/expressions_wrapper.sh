#!/bin/bash
#SBATCH -p batch              # Standard CPU partition
#SBATCH -N 1                  # 1 Node
#SBATCH -n 1                  # 1 Task
#SBATCH -c 1                  # 1
#SBATCH --mem=2G
#SBATCH --time=00:05:30       # Adjust this based on your chain length test
#SBATCH -J run_8       # Job name
#SBATCH -o /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/logs/run_8       # Standard output log

if [ $# -ne 2 ]; then
  echo "Usage: $0 <run_number> [base_root]"
  exit 1
fi

python3 /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/scripts/run_expressions.py --run $1 $2