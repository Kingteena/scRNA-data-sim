#!/bin/bash
#SBATCH -p batch              # Standard CPU partition
#SBATCH -N 1                  # 1 Node
#SBATCH -n 1                  # 1 Task
#SBATCH -c 1                  # 1
#SBATCH --mem=400M
#SBATCH --time=02:40:00       # Adjust this based on your chain length test
#SBATCH -J beast        # Job name
#SBATCH -o /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast/log        # Standard output log


if [ $# -ne 1 ]; then
  echo "Usage: $0 <XML FILE>"
  exit 1
fi

module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
module load beagle-lib/4.0.0-GCC-11.2.0

beast -beagle_CPU -beagle_SSE $1
