#pip install Biopython

import os
import re
import numpy as np
import pandas as pd
from Bio import Phylo

def get_cellnames(path):
    """
    Reads Newick tree from Cellcoal tree.000X file and gets the name of the cells
    """
    tree = Phylo.read(path, "newick")
    leaves = tree.get_terminals()
    names = [leaf.name for leaf in leaves[:-1]]
    return names

def assign_cells_to_states(cnames, nstates):
    """
    Assigns cells to one of k states using a Dirichlet distribution.

    Parameters:
    - cnames (list): contains cell names for a tree.
    - alpha (list or array): Dirichlet concentration parameters for the states.

    Returns:
    - cell_states: A list containing the State of each cell.
    """
    alpha_states = np.ones(nstates) # Symmetric dirichlet; can be adjusted

    # Make sure there is at least 2 distinct cell assignments
    while True:
        # Generate Dirichlet distribution for probabilities across states
        state_probabilities = np.random.dirichlet(alpha_states, size=len(cnames))
            
        # Assign each cell to a state based on the highest probability
        cell_assignments = np.argmax(state_probabilities, axis=1)

        if len(set(cell_assignments)) >= 2:
            break
    
    return cell_assignments

def simulate_state_ges(nstates, ngenes, alpha_ges):
    """
    Simulate Gene Expression States (GES) for each state.
    
    Parameters:
        nstates: int, number of cell states (k)
        ngenes: int, number of genes (m)
        alpha_ges: float, concentration parameter for the Dirichlet distribution
        
    Returns:
        ges: array of shape (nstates, ngenes), GES for each state
    """
    ges = np.random.dirichlet([alpha_ges] * ngenes, nstates)
    return ges

def simulate_cell_counts(cell_assignments, ges, minlib, maxlib):
    """
    Simulate counts for cells based on their assigned state.
    
    Parameters:
        cell_assignments: array of shape (ncells,), state assignments for each cell
        ges: array of shape (num_states, ngenes), GES for each state
        minlib (int): minimum library size
        maxlib (int): maximum library size
        
    Returns:
        expression_matrix: data_frame (ncells, ngenes) with simulated counts
    """
    ncells = len(cell_assignments)
    ngenes = ges.shape[1]
    library_sizes = np.random.uniform(minlib, maxlib, ncells)  # Random library sizes
    expression_matrix = np.zeros((ncells, ngenes), dtype=int)
    
    for i, state in enumerate(cell_assignments):
        state = ges[state]                      # Get GES for the cell's state
        scaling_factor = library_sizes[i]
        gene_means = state * scaling_factor
        expression_matrix[i, :] = np.random.poisson(gene_means)

    #create data_frame
    genes = [f"Gene{i+1}" for i in range(ngenes)]                        # Gene names
    cells = [f"tumcell{str(i+1).zfill(4)}" for i in range(ncells)]  # Cell names
    expression_matrix = pd.DataFrame(expression_matrix, index=cells, columns=genes)
    
    return expression_matrix

def parse_vcf(vcf_file):
    """Parse the VCF file into a DataFrame."""
    with open(vcf_file, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('##')]
    header = lines[0].split('\t')
    data = [line.split('\t') for line in lines[1:]]
    vcf_df = pd.DataFrame(data, columns=header)
    #del vcf_df[vcf_df.columns[-1]]  delete last column if needed
    return vcf_df

def update_vcf_with_expression(vcf_df, counts_df, genome_length):
    """
    Update the VCF file based on expression data.

    Parameters:
    - vcf_df: DataFrame containing the VCF data.
    - counts_df: DataFrame containing gene expression data (rows: cells, columns: genes).
    - genome_length: Total length of the genome, used to map positions to genes.

    Returns:
    - Updated VCF DataFrame.
    """
    ngenes = counts_df.shape[1]            # Number of genes
    gene_length = genome_length // ngenes  # Length of each gene

    for idx, row in vcf_df.iterrows():
        # Determine which gene this SNV belongs to based on its position
        position = int(row['POS'])
        gene_index = position // gene_length
        
        if gene_index >= ngenes:
            continue  # Skip if position is out of range

        # Get the gene name from counts_df columns
        gene = counts_df.columns[gene_index]

        for cell in counts_df.index:            # Iterate over cells
            if counts_df.loc[cell, gene] == 0:  # Check if expression is zero
                # Update the genotype for this cell in the VCF
                genotype = row[cell]
                updated_genotype = re.sub(r'\b[01]\|[01]:', '.|.:', genotype)
                vcf_df.at[idx, cell] = updated_genotype
    
    return vcf_df


def write_vcf(vcf_df, output_vcf):
    """Write the updated VCF DataFrame to a file.
    
    Parameters:
    - vcf_df: DataFrame containing the updated VCF data.
    - output_vcf: Path to save the updated VCF file.
    """
    with open(output_vcf, 'w') as f:
        f.write("##fileformat=VCFv4.3\n")
        f.write("##source=Updated with zero-expression filtering\n")
        f.write('\t'.join(vcf_df.columns) + '\n')
        for _, row in vcf_df.iterrows():
            f.write('\t'.join(row) + '\n')


#if __name__ == "__main__":
#    args = parse_args()
#    vcf2xml(args.vcf, args.xml_template, args.output_xml, args.encoding)

###  Parameters  ###
# genome_length = length of the dna Sequence
# ngenes = number of genes
# nstates = number of states 
# alpha_states = np.ones(nstates) # Symmetric dirichlet; can be adjusted
# alpha_ges = 0.5                 # Dirichlet concentration parameter
# maxlib = 60000
# minlib = 40000


### file paths ###
# vcf_file = "results/vcf_dir/vcf.0001"
# updated_vcf_file = "results/vcf_dir/scRNAvcf.0001"


