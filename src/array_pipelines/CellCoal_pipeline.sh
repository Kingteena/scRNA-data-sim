#!/bin/bash
#SBATCH --job-name=cellcoal
#SBATCH --output=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output/cellcoal/cellcoal_logs/run_%a.log    # Saves console output to cellcoal_logs/run_1.log, etc.
#SBATCH --array=1-8                # Launches 8 independent tasks
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1                       # CellCoal is single-threaded
#SBATCH --mem=210GB                # Adjust based on your genome size. sequencing errors and number of cells 
#SBATCH --time=01:30:00            # Adjust based on expected runtime

# Notification configuration
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahana.meenachisundaram@student.adelaide.edu.au

# CellCoal Simulation Script 
# The script loops through an set array, which indicates the number of different parameter files that need to be used to simulate cancer cell lineages 
# Within this project, 8 different combination of parameters were needed, thus the "--array-1-8", lauching 8 independant tasks 
# The script runs cellcoal into the specified output folders, as specified within the script and parameter file
# sbatch Cellcoal_SimScript to launch the script, after making the adjustments highlighted by the comments throughout the script 

# Setting Up Paths 

PROJROOT=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation  #Set your own Project Root
CELLCOAL_BIN=${PROJROOT}/tools/cellcoal-1.1.1 
PARAMETERS_SOURCE="parameters/run_${SLURM_ARRAY_TASK_ID}.txt" #Setting location to the parameter files 
OUTPUT=${PROJROOT}/output
CELLCOAL_FOLDER=${OUTPUT}/cellcoal
RUN_DIR="run_${SLURM_ARRAY_TASK_ID}"
CELLCOAL_TEMPLATE=/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/scripts/templates/cellcoal_template.sh  # Adjust for your system, use direct path

# Preparing the enviroment for the run interation 
mkdir -p $CELLCOAL_FOLDER
mkdir -p "$CELLCOAL_FOLDER/$RUN_DIR"
mkdir -p $CELLCOAL_FOLDER/cellcoal_logs 

# Move into the run directory so each run stays isolated
cd "$CELLCOAL_FOLDER/$RUN_DIR" 

echo "Starting simulation for task ${SLURM_ARRAY_TASK_ID}..."

# Running Cellcoal: 
$CELLCOAL/bin/cellcoal-1.1.0 "-F$PARAMETERS_SOURCE"

echo "Task ${SLURM_ARRAY_TASK_ID} complete."

echo "Finished with Simulation of Cancer Cell Lineages" 
echo "----------------------------------------------------------------------------------------------------------------------------------"

