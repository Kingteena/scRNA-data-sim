# Cancer Cell Lineage Simulation Project
---
This project aims to investigate the extent to which sequencing noise impacts the accuracy of the reconstructed phylogeny using Maximum Likelihood (CellPhy) and Bayesian Models (BEAST2)

To Get started, clone the GitHub Repo and follow the steps in [Usage & Installation](https://github.com/Kingteena/scRNA-data-sim/wiki/Installation-%26-Setup)

For usage, see [scripts](https://github.com/Kingteena/scRNA-data-sim/wiki/Scripts)


For an explanation about the theory, read [this](https://github.com/Kingteena/scRNA-data-sim/wiki/Ground-Truth-Generation)


## Tools Used
- [CellCoal](https://github.com/dapogon/cellcoal/)  
- [CellPhy](https://github.com/amkozlov/cellphy/)
- [BEAST2](http://www.beast2.org/)
    - [beast-phylonco](https://github.com/bioDS/beast-phylonco) - Available as a package in BEAST2
    - [BEAGLE](https://github.com/beagle-dev/beagle-lib) - optional but recommended for better performance
- [Conda](https://docs.conda.io/en/latest/) - Optionally use Mamba on HPC systems for faster package resolution
- [TreeDist](https://github.com/ms609/TreeDist) - R Package to compare similarity between phylogenetic trees


## Pipeline
The simulation pipeline is organized into the following steps:
1. **Simulation of Cancer Cell Lineages**: Simulate the cancer cell lineages via using CellCoal
2. **Data Formatting**: Convert CellCoal output into formats compatible with downstream analysis tools
3. **Phylogenetic Inference**: Use CellPhy and BEAST2 with the phylonco package to infer phylogenetic trees from the simulated data.
4. **Analysis and Visualization**: Analyze the inferred trees and visualize the results 

For more details on the scripts used in each step, see the [Wiki](https://github.com/Kingteena/scRNA-data-sim/wiki)