#!/bin/bash

TEMPLATE_SCRIPT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/beast_template.sh
OUTPUT_DIR=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/beast

# Take the VCF file as an argument
if [ $# -ne 1 ]; then
  echo "Usage: $0 <VCF FILE>"
  exit 1
fi

VCF_FILE=$1
FILE_NAME=$(basename ${VCF_FILE%.vcf})
mkdir -p ${OUTPUT_DIR}/$FILE_NAME/xml
mkdir -p ${OUTPUT_DIR}/$FILE_NAME/run
mkdir -p ${OUTPUT_DIR}/$FILE_NAME/tree

# Create 2 XML files for BEAST analysis
echo "Creating XML files for BEAST analysis from $FILE_NAME"

python3 /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/vcf2xml.py $VCF_FILE /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/test_gt16_error.xml ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml --filtering_density 0.3

# Run BEAST on both XML files
echo "Running BEAST on the generated XML files" 

sed 's|$(seed)|'${OUTPUT_DIR}'/'${FILE_NAME}'/run/1|g' ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml > ${OUTPUT_DIR}/${FILE_NAME}/xml/1.xml
sed 's|$(seed)|'${OUTPUT_DIR}'/'${FILE_NAME}'/run/2|g' ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml > ${OUTPUT_DIR}/${FILE_NAME}/xml/2.xml
sed 's|$(seed)|'${OUTPUT_DIR}'/'${FILE_NAME}'/run/3|g' ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml > ${OUTPUT_DIR}/${FILE_NAME}/xml/3.xml
sed 's|$(seed)|'${OUTPUT_DIR}'/'${FILE_NAME}'/run/4|g' ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml > ${OUTPUT_DIR}/${FILE_NAME}/xml/4.xml


rm ${OUTPUT_DIR}/${FILE_NAME}/xml/template.xml

sbatch $TEMPLATE_SCRIPT ${OUTPUT_DIR}/${FILE_NAME}/xml/1.xml
sbatch $TEMPLATE_SCRIPT ${OUTPUT_DIR}/${FILE_NAME}/xml/2.xml
sbatch $TEMPLATE_SCRIPT ${OUTPUT_DIR}/${FILE_NAME}/xml/3.xml
sbatch $TEMPLATE_SCRIPT ${OUTPUT_DIR}/${FILE_NAME}/xml/4.xml

echo "BEAST analysis jobs submitted for $FILE_NAME"

# Take trees files and summarize them using TreeAnnotator
# echo "Summarizing BEAST trees using TreeAnnotator"

# module load Beast/2.7.7-GCC-12.3.0-CUDA-12.4.1
# treeannotator -topology CCD0 ${OUTPUT_DIR}/${FILE_NAME}/run/1.trees ${OUTPUT_DIR}/${FILE_NAME}/tree/1.tree
# treeannotator -topology CCD0 ${OUTPUT_DIR}/${FILE_NAME}/run/2.trees ${OUTPUT_DIR}/${FILE_NAME}/tree/2.tree 
