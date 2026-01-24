# # Paths to files
# TREE_INPUT_PATH="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellcoal_simulation_outputs/u9_c10_E5_D6/trees_dir/trees.0001"
# VCF_INPUT_PATH="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellcoal_simulation_outputs/u9_c10_E5_D6/vcf_dir/0001.vcf"               # Path to the input VCF file
# VCF_OUTPUT_PATH="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/cellcoal_simulation_outputs/u9_c10_E5_D6/RNAvcf.0001"            # Path to save the updated VCF file
# EXPRESSION_PROFILES_PATH="expression_profiles.csv"               # Path to save expression profiles

# TODO: CONVERT THIS SCRIPT TO PYTHON SCRIPT FOR BETTER READABILITY AND MAINTAINABILITY

# Take three arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <tree_input_path> <vcf_input_path> <vcf_output_path>"
    exit 1
fi

TREE_INPUT_PATH="$1"
VCF_INPUT_PATH="$2"
VCF_OUTPUT_PATH="$3"

# Parameters
GENOME_LENGTH=500000                    # Total genome length
NUM_STATES=5                            # Number of cell states
NUM_CELLS=10                            # Number of cells for expression profiles
NUM_GENES=30000                        # Number of genes for expression profiles
ALPHA_GES=0.5                           # Dirichlet concentration parameter 
MINLIB=40000                            # Minimum library size
MAXLIB=60000                            # Maximum library size

# Create expression profiles
echo "Generating expression profiles and updating vcf file..."
python3 -c "
cnames = get_cellnames(path='${TREE_INPUT_PATH}')
cell_assignments = assign_cells_to_states(cnames, nstates=${NUM_STATES})
ges = simulate_state_ges(nstates=${NUM_STATES}, ngenes=${NUM_GENES}, alpha_ges=${ALPHA_GES})
counts_df = simulate_cell_counts(cell_assignments, ges, minlib=${MINLIB}, maxlib=${MAXLIB})
vcf_df = parse_vcf(vcf_file='${VCF_INPUT_PATH}')
vcf_df = update_vcf_with_expression(vcf_df, counts_df, genome_length=${GENOME_LENGTH})
write_vcf(vcf_df, output_vcf='${VCF_OUTPUT_PATH}')
"

# convert vcf to xml for beast 
# echo "Converting vcf to xml..."
# python3 python/vcf2xml.py data/vcf_exp/RNAvcf.0001 data/xml/cellcoal_temp.xml data/xml/cellcoal_0001.xml

echo "VCF updates with expressions successfully!"
