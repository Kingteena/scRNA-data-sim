#!/bin/bash
#SBATCH -p batch              # Standard CPU partition
#SBATCH -N 1                  # 1 Node
#SBATCH -n 1                  # 1 Task
#SBATCH -c 1               # 4
#SBATCH --mem=4G              # 4GB Total RAM (2GB for Java + 2GB buffer for BEAGLE/OS)
#SBATCH --time=00:20:00       # Adjust this based on your chain length test
#SBATCH -J beast        # Job name
#SBATCH -o instances.log          # Standard output log


if [$# -ne 1]; then
  echo "Usage: $0 <XML FILE>"
  exit 1
fi

module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
module load beagle-lib/4.0.0-GCC-11.2.0

beast -beagle_CPU -instances 2 $1
