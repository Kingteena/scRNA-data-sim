#!/bin/bash
# Creates a summary table of folder name vs SNV count from all the replicates generated in the simulation 
# Also generates an average SNV count from all the replicates generated 

#please change the directory as per your system and parameters (replicates value)
SUMMARY_FILE=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/Simulated_SNV_summary.txt
TOTAL_SNVs=0
FILE_COUNT=0
REPLICATES=30

touch $SUMMARY_FILE
echo -e "Summary Table of the SNVs generated within the Simulated Replicates\n" >> $SUMMARY_FILE
echo -e "VCF_File\t|\tSNV_Count" >> $SUMMARY_FILE
echo "----------------------------------------" >> $SUMMARY_FILE

for folder in /hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/; do
    echo "Processing folder: $(basename $folder)" 
    for file in $folder/cleaned_reps/* ; do
        SNV=$(grep -v "^#" "$file" | wc -l)
         
         # Add the count to the total value
         TOTAL_SNVs=$((TOTAL_SNVs + SNV))

        # Increment the file count
        FILE_COUNT=$((FILE_COUNT + 1))

        echo -e "$(basename $file)\t|\t$SNV" >> $SUMMARY_FILE
    done
done

# Calculate the average SNV value (use 'bc' for floating point division)
            if [[ "$FILE_COUNT" == "$REPLICATES" ]]; then
            AVERAGE_SNVs=$(echo "scale=2; $TOTAL_SNVs / $FILE_COUNT" | bc)
            echo -e "\nThe Average number of SNVs generated per simulated replicate is $AVERAGE_SNVs" >> $SUMMARY_FILE
            else
            AVERAGE_SNVs=0
            fi

echo "Summary Table of the Filtered simulated Cancer Cell Lineages has been generated in the cellcoal folder of outputs"