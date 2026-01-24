#!/bin/bash
#SBATCH -p batch              # Standard CPU partition
#SBATCH -N 1                  # 1 Node
#SBATCH -n 1                  # 1 Task
#SBATCH -c 1                  # 4
#SBATCH --mem=4G              # 4GB Total RAM (2GB for Java + 2GB buffer for BEAGLE/OS)
#SBATCH --time=10:00:00       # Adjust this based on your chain length test
#SBATCH -J beast        # Job name
#SBATCH -o logs/%x_%j.log          # Standard output log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a1930250@adelaide.edu.au

module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
module load beagle-lib/4.0.0-GCC-11.2.0i

beast -beagle_CPU -beagle_SSE  /hpcfs/users/a1930250/beast/<XML FILE>
