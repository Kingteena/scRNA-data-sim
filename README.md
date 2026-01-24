# Cancer Cell Lineage Simulation Project

**TODO: Insert overview of the data simulation project here.**

**NOTE:** This Project was designed to be run on the [Phoenix HPC cluster](https://www.adelaide.edu.au/technology/research/high-performance-computing/phoenix-on-premise-hpc) at the University of Adelaide. Details have been given on how to run on other systems, however adjustments may be necessary. 

## Table of Contents
- [Tools Used](#tools-used)
- [Installation](#installation)
- [Pipeline](#pipeline)
- [Scripts](#scripts)


## Tools Used
- [CellCoal](https://github.com/dapogon/cellcoal/)  
- [CellPhy](https://github.com/amkozlov/cellphy/)
- [BEAST2](http://www.beast2.org/)
    - [beast-phylonco](https://github.com/bioDS/beast-phylonco) - Available as a package in BEAST2
    - [BEAGLE](https://github.com/beagle-dev/beagle-lib) - optional but recommended for better performance
- [Conda](https://docs.conda.io/en/latest/) - Optionally use Mamba on HPC systems for faster package resolution
- [TreeDist](https://github.com/ms609/TreeDist) - R Package to compare similarity between phylogenetic trees

## Installation

All the required tools should be installed into the `tools` directory.

Note: The installation instructions below are tailored for HPC systems using `module` to manage software. On other systems, simply install from the respective sources or use package managers like `conda` or `brew`.

### Conda Environment
Conda is recommended for managing the software dependencies. Conda can be directly loaded on HPC systems using `module`:

```bash
module load Anaconda3
```
You can then update conda, create a conda environment using the provided `environment.yml` file and update the environment as follows:

```bash
conda update conda
conda env create -f environment.yml
conda activate Simulation
conda update --all
```

### CellCoal
You can then install CellCoal by following the instructions on their respective GitHub repositories. For this Simulation, we employed the version CellCoal-1.1.1. 
 
Cellcoal needs to be complied after downloading:

```bash
cd cellcoal-x.y.z
make
```

### CellPhy
CellPhy can be used directly from the provided script, in this simulation project, cellphy-0.9.4 was employed. However, installing it this way will not allow it to generate .svg visualizations of the trees. This requires installing R and the required packages as described in the CellPhy repository. However, the R packages are already included in the conda environment. As such, R needs to be installed separately if you wish to generate the visualizations.

To install R, you can use `module` or `conda`. Using `module` on HPC systems is recommended:

```bash
# On HPC systems:
module load R

# Or using conda:
conda install R
```

### BEAST2 and phylonco package
When installing Beast2 it is also recommended to install [beagle](https://github.com/beagle-dev/beagle-lib) separately for better performance with BEAST2. Both can be loaded via `module` on HPC systems.

```bash
module load Beast beagle-lib
```

To install the phylonco package for BEAST2, you can use the BEAST2 package manager:
(On local systems, the BEAUTi GUI can be used to install packages as well.)

```bash
packagemanager -add phylonco
```
### TreeDist R Package 
This is an R package developed by Professor Martin R. Smith of Durham University, that utilises various tools to quantify the topological distances between unweighted phylogenetic trees. In this project it was employed to calculate the Generalised Robinson-Foulds distance between two or more trees, to reflect the similarity between them. 

Install and load the library using a CRAN library with the below commands
```R
install.packages('TreeDist')
library('TreeDist')
```
You can install the development version of the package with:
```R
if(!require("curl")) install.packages("curl")
if(!require("remotes")) install.packages("remotes")
remotes::install_github("ms609/TreeDist")
```

## Pipeline
The simulation pipeline is organized into the following steps:
1. **Simulation of Cancer Cell Lineages**: 

- The cancer cell lineages were simulated using **CellCoal (v1.1.1)**, with a loose control of ~52000 SNVs per replicate generated, as this closely replicates the value of SNVs present in practically found scRNA-seq data. 
- The exact parameters selected for the simulation, were acheived through trial and error and are found within a parameter file **parameters_cellcoal.txt** 
- Within the key control parameters implemented, the SNVs were introduced predominately by somatic mutation rate (u), germline SNP rate (c), fixed allelic dropout (D) and sequencing error (E). 

2. **Data Formatting**: 


Convert CellCoal output into formats compatible with downstream analysis tools

3. **Phylogenetic Inference**: Use CellPhy and BEAST2 with the phylonco package to infer phylogenetic trees from the simulated data.
4. **Analysis and Visualization**: Analyze the inferred trees and visualize the results 

## Scripts
For more details about thes scripts, please refer to the respective [README](scripts/README.md).
For more detailed instructions and examples, please refer to the individual scripts and configuration files included in this repository.

