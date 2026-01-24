#!/bin/bash
# This script generate the parameter file for the specified ranges and calls the slurm script to run the simulations
# Parameters:
# Somatic mutation rate: 2e-6 to 2e-10 (u = 6 to 10)
# Germline mutation rate: 2e-6 to 2e-10 (c = 6 to 10)
# Sequencing error rate: 1e-3 to 1e-5 (E = 3 to 5)
# Allelic dropout rate: 1e-3 to 1e-6 (D = 3 to 6)

# It generates a file with the name parameters_u{u}_c{c}_E{E}_D{D}.txt
# where u, c, E, D are the respective parameter exponents

# We want c<u, D<E
# as we're dealing with negative integers, we want $c > $u and $D > $E

TEMPLATE_PARAMETER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellcoal-1.1.1/parameters_template.txt
TEMPLATE_SCRIPT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/cellcoal_template.sh
BASEDIR=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation
OUTPUT_FOLDER=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal

for u in {1,5,10}; do
    for c in {1,5,10}; do
        if ((c > u)); then
            for E in {1,5,10}; do
                D=5
                    # if ((D > E )); then # We need D < E but as its

                        PREFIX="u${u}_c${c}_E${E}_D${D}"
                        

                        mkdir -p ${OUTPUT_FOLDER}/${PREFIX}
                        
                        echo "Generating parameter file: ${PREFIX}"

                        # Create the template parameter file with the desired parameters
                        sed -e "s/<Som_rate>/${u}e-9/g" \
                            -e "s/<Germ_rate>/${c}e-10/g" \
                            -e "s/<Seq_error>/${E}e-4/g" \
                            -e "s/<ADO_rate>/1e-${D}/g" \
                            -e "s|<Folder_name>|$OUTPUT_FOLDER/${PREFIX}|g" \
                            $TEMPLATE_PARAMETER > $OUTPUT_FOLDER/${PREFIX}/parameters.txt

                        # $TEMPLATE_SCRIPT $OUTPUT_FOLDER/${PREFIX}/parameters.txt $OUTPUT_FOLDER/${PREFIX}
                        sbatch $TEMPLATE_SCRIPT $OUTPUT_FOLDER/${PREFIX}/parameters.txt $OUTPUT_FOLDER/${PREFIX}
                        sleep 5 # to avoid overloading the scheduler and being rate limited
                    # fi
                # done
            done
        fi
    done
done


