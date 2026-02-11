#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --time=00:30:00
#SBATCH --mem=72GB
#SBATCH --output="logs/%x-%j.out"
#SBATCH --error="logs/%x-%j.err"

if [ $# -ne 2 ]; then
  echo "Usage: $0 <parameter_file> <output_folder>"
  exit 1
fi

#Key areas of directory present and creating them within the folder 
PROJROOT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation
CELLCOAL=${PROJROOT}/tools/cellcoal-1.1.1
PARAMETERS=$1
OUTPUT=$2
Folder_name=output

#Running Cellcoal and Simulating the Cancer Cell Lineages 
echo "Now running Cellcoal and Simulating Cancer cell lineages"

${CELLCOAL}/bin/cellcoal-1.1.0 -F${PARAMETERS} 

echo "Done with running CellCoal"
