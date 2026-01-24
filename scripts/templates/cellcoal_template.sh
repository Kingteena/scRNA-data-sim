#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --time=00:30:00
#SBATCH --mem=72GB

#SBATCH --output="logs/%x-%j.out"
#SBATCH --error="logs/%x-%j.err"

module load BCFtools/1.17-GCC-11.2.0

if [ $# -ne 2 ]; then
  echo "Usage: $0 <parameter_file> <output_folder>"
  exit 1
fi

#Key areas of directory present
PROJROOT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation
CELLCOAL=${PROJROOT}/tool/cellcoal-1.1.1
OUTPUT=$2
FILTERED=$OUTPUT/cleaned_reps
PARAMETERS=$1


## Check the project root exists
if [[ -d ${PROJROOT} ]]; then
  echo -e "Found ${PROJROOT}\n"
else
  echo -e "${PROJROOT} not found.\nExiting with Error code 1"
  exit 1
fi

mkdir -p $OUTPUT

#Running Cellcoal and Simulating the Cancer Cell Lineages 
echo "Now running Cellcoal and Simulating Cancer cell lineages"

${CELLCOAL}/bin/cellcoal-1.1.0 -F${PARAMETERS}

echo "Done with CellCoal"

mkdir -p $FILTERED

# cd ${CELLCOAL}/$Folder_name/vcf_dir  

# renaming the file to have .vcf suffix 
# ##*. remove everything after the pattern 
  for f in $OUTPUT/vcf_dir/vcf.*; do
    mv "$f" "$OUTPUT/vcf_dir/${f##*.}.vcf"
  done 

echo "Now filtering the files"

# remove the healthy cell or outgroup and forms a filtered file version
  for file in $OUTPUT/vcf_dir/*.vcf; do
    bcftools view -s ^healthycell "$file" -o "${FILTERED}/$(basename $file)"
done

echo "Finished with Simulation of Cancer Cell Lineages" 
