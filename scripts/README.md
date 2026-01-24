# Data Simulation Scripts

This directory contains scripts used for simulating data for various analyses. This file provides an overview of the available scripts and their purposes.

# Templates

These scripts serve as templates for running simulations with different parameters. Users can modify these templates to suit their specific needs via command-line arguments. They are also written in a way that can be called with `sbatch <script_name>` for parallel execution on a computing cluster. They can also simply be called directly from the command line with `./<script_name>` if desired. 

These are usually called from other scripts that set up the parameters for multiple runs.

These scripts include:
- `cellcoal_template.sh`:  
    - A template for running CellCoal simulations with specified parameters. It sets up the environment and calls the CellCoal tool with the provided parameters.
    - This script changes somatic and germline mutation rate, sequencing error rate and allelic dropout rate. The rest of the parameters are fixed within the parameters file in cellcoal_parameters.txt.
    - It then takes the output from CellCoal and removes 'healthycell' samples to remove the outgroup. This is then stored in a separate folder called cleaned_reps within the output directory.
    - Usage: `sbatch cellcoal_template.sh <parameter_file> <output_folder>`
- `cellphy_template.sh`:
    - A template for running CellPhy analyses on VCF files generated from CellCoal simulations.
    - This script takes VCF files from the cleaned_reps directory and runs CellPhy to infer phylogenetic trees.
    - Usage: `sbatch cellphy_template.sh <input_vcf_file> <output_directory>`

## Logging
Both template scripts log their output to files named `scriptname_timestamp.log` to `logs/` directory within the output directory specified. This helps in tracking the progress and debugging any issues that may arise during execution.

## Changing Parameters on HPC
When changing parameters, ensure that the sbatch options (e.g., memory, time) are appropriate for the new settings. 

We have noticed that CellCoal is memory hungry, and for 50k SNVs, 5 replicates take about 20GB of RAM, and it scales logarithmically with the number of replicates at about 70 GB for 50 replicates.

Meanwhile, CellPhy is less memory hungry (often <1GB) but more time consuming. However, it is embarrassingly parallel, so you can run many replicates in parallel to speed things up.

# Callers
These scripts are designed to call the template scripts with specific parameters. They often loop over a set of conditions or configurations to automate the simulation process across multiple scenarios. These scripts are typically executed with `./<caller name>.`

**NOTE:**: These scripts are desined to run on an HPC cluster by calling template scripts via `sbatch`. To run them on a local machine, you may need to modify the calls to the template scripts to execute them directly (i.e., `./<template_script>` instead of `sbatch <template_script>`).