#!/bin/bash

# Path to your binary
CELLPHY_BIN="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/tools/cellphy-0.9.4/bin/raxml-ng-cellphy-linux"

echo "scanning .vcf files..."
printf "%-30s | %-10s | %-10s\n" "FILENAME" "THREADS" "MEMORY"
echo "-------------------------------------------------------------"

# Use nullglob to avoid issues if no .vcf files are present
shopt -s nullglob
for vcf in *.vcf; do
    # Run the parse command, capture output, suppress the binary file creation if possible
    # We grep strictly for the lines we need to avoid spam
    OUTPUT=$($CELLPHY_BIN --parse --msa "$vcf" --model GT16+FO --threads 1 2>&1)

    # Extract the numbers using awk/grep magic
    THREADS=$(echo "$OUTPUT" | grep "Recommended number of threads" | awk -F': ' '{print $2}')
    MEM=$(echo "$OUTPUT" | grep "Estimated memory requirements" | awk -F': ' '{print $2}')

    # Print nicely
    printf "%-30s | %-10s | %-10s\n" "$vcf" "$THREADS" "$MEM"
    
    # Cleanup: --parse creates a binary .rba file you don't need yet. Delete it to save space.
    rm -f "${vcf}.raxml.rba"
done