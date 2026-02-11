"""
Expression Filter for VCF Files.

This module integrates gene expression simulation with VCF genotype filtering.
It simulates gene expression counts based on cell states and masks variants
in the VCF that fall within genes with zero expression.
"""
import os

import numpy as np
import pandas as pd


def get_cellnames(vcf_df: pd.DataFrame) -> list[str]:
    """
    Reads headers of the vcf dataframe to return a list of the name of the cells
    """
    return vcf_df.columns[9:].tolist()


def assign_cells_to_states(cnames: list[str], nstates: int) -> np.ndarray[int]:
    """
    Assigns cells to one of k states using a Uniform distribution. 
    Ensures at least 2 distinct states are chosen, unless nstates = 1.

    Parameters:
    - cnames: contains cell names for a tree.
    - nstates: Number of states (k) to assign cells to.

    Returns:
    - cell_states: A list containing the State of each cell.
    """
    while True:
        cell_assignments = np.random.randint(0, nstates, size=len(cnames))

        if nstates == 1:
            print("Warning: Only one state specified. All cells will be assigned to the same state.")
            break
        # Ensure at least 2 distinct states are chosen
        if len(set(cell_assignments)) >= 2:
            break

    return cell_assignments


def simulate_state_ges(
    nstates: int, ngenes: int, alpha_ges: list[float]
) -> np.ndarray[np.ndarray[float]]:
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


def simulate_cell_counts(
    cell_assignments: np.ndarray[int],
    ges: np.ndarray[np.ndarray[float]],
    cnames: list[str],
    minlib: int,
    maxlib: int,
) -> pd.DataFrame:
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

    # Generate Library Sizes (Vectorized)
    library_sizes = np.random.uniform(minlib, maxlib, ncells)

    # Expand GES to (ncells, ngenes)
    cell_gene_means = ges[cell_assignments]

    # Reshape library_sizes to (ncells, 1) to broadcast across genes
    means = cell_gene_means * library_sizes[:, None]
    expression_matrix = np.random.poisson(means)
    genes = [f"Gene{i+1}" for i in range(ngenes)]

    return pd.DataFrame(expression_matrix, index=cnames, columns=genes)


def parse_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Parses VCF file of the given path into a DataFrame, filtering out comment lines and outgroup
    and ensuring POS is integer.
    """
    with open(vcf_file, "r", encoding="utf-8") as f:
        # Single pass: Filter comments, strip, and split all at once
        data = [line.strip().split("\t") for line in f if not line.startswith("##")]

    df = pd.DataFrame(data[1:], columns=data[0])
    df["POS"] = df["POS"].astype(int)  # Ensure POS is integer for later processing

    # Delete the 'healthycell' column if it exists, as we want to remove the outgroup
    if "healthycell" in df.columns:
        df.drop(columns=["healthycell"], inplace=True)
    return df


def update_vcf_with_expression(
    vcf_df: pd.DataFrame, counts_df: pd.DataFrame, genome_length: int
) -> pd.DataFrame:
    """
    Optimized update using Boolean Masking.
    Includes safety check for POS column type.
    """
    ngenes = counts_df.shape[1]

    if ngenes == 0:
        raise ValueError("Number of genes cannot be zero.")
    elif ngenes > genome_length:
        raise ValueError(
            f"Number of genes ({ngenes}) cannot be greater than genome length ({genome_length})."
        )

    gene_length = genome_length // ngenes

    # --- FIX IS HERE ---
    # Force POS to be integers. If it's already int, this does nothing (fast).
    # If it's string, it fixes it instantly.
    # vcf_df['POS'] = vcf_df['POS'].astype(int)
    # -------------------

    # 1. Pre-calculate which gene every variant belongs to
    variant_gene_indices = vcf_df["POS"].values // gene_length

    # 2. Loop ONLY over cells (columns)
    valid_cells = [c for c in counts_df.index if c in vcf_df.columns]

    for cell in valid_cells:
        zero_gene_indices = np.where(counts_df.loc[cell].values == 0)[0]

        if len(zero_gene_indices) == 0:
            continue

        mask = np.isin(variant_gene_indices, zero_gene_indices)

        if np.any(mask):
            vcf_df.loc[mask, cell] = vcf_df.loc[mask, cell].str.replace(
                r"^[01]\|[01]:", ".|.:", regex=True
            )

    return vcf_df


def write_vcf(vcf_df: pd.DataFrame, output_vcf: str) -> None:
    """Write the updated VCF DataFrame to a file.

    Parameters:
    - vcf_df: DataFrame containing the updated VCF data.
    - output_vcf: Path to save the updated VCF file.
    """
    vcf_df["POS"] = vcf_df["POS"].astype(str)  # Ensure POS is string for output

    with open(output_vcf, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.3\n")
        f.write("##source=Updated with zero-expression filtering\n")
        f.write("\t".join(vcf_df.columns) + "\n")

        for row in vcf_df.itertuples(index=False, name=None):
            f.write("\t".join(row) + "\n")


def process_single_file(vcf_input_path: str, output_path: str, config: dict) -> str:
    """
    The worker task that handles ONE file completely.
    """
    try:
        # Input vcf and parse it
        vcf_df = parse_vcf(vcf_input_path)
        cnames = get_cellnames(vcf_df=vcf_df)

        # Calculate the variables
        cell_assignments = assign_cells_to_states(cnames, nstates=config["NUM_STATES"])
        ges = simulate_state_ges(
            config["NUM_STATES"], config["NUM_GENES"], config["ALPHA_GES"]
        )
        counts_df = simulate_cell_counts(
            cell_assignments, ges, cnames, config["MINLIB"], config["MAXLIB"]
        )

        # Update VCF
        vcf_df = update_vcf_with_expression(
            vcf_df, counts_df, genome_length=config["GENOME_LENGTH"]
        )

        # Write updated VCF
        write_vcf(vcf_df, output_vcf=output_path)

        return f"SUCCESS: {os.path.basename(vcf_input_path)}"

    except Exception as e:
        return f"Error {os.path.basename(vcf_input_path)}: {str(e)}"
