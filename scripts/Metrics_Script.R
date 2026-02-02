
## Final Step of Simulation Workflow
## R script summarises the accuracy of phylogentic tree inferences from CellPhy and BEAST through key comparison  metrics
## BoxPlot (.pdf) and Table (.csv) are generated for each replicates set of three trees inputted from script caller
## FOR ML, CellPhy outputs unrooted trees, for its Best.Tree file, thus the true tree or Cellcoal tree is also unrooted in metric calcs 
## FOR BAYESIAN, BEAST already outputs rooted trees for its tree file


# trailingOnly = TRUE: ignores all the background "system" info and only focuses on specific arguments you typed in Bash in the caller 
# args:  Creates a vector of each argument as a numbered item 

args <- commandArgs(trailingOnly = TRUE)

TRUETREE_PATH <- args[1]
BEASTtree_PATH <- args[2]
CELLPHYtree_PATH <- args[3]
OUTPUT_PATH <- args[4]

#TRUETREE_PATH="/home/sahana/output/cellcoal/tree_files/0003.trees"
#BEASTtree_PATH = "/home/sahana/output/beast/0003/tree/0003.tree"
#CELLPHYtree_PATH="/home/sahana/output/cellphy/0003/0003.vcf.raxml.bestTree"

message("Loading in any required libraries into current workspace")
# Loading in any required libraries into current workspace 

if (!require("ape")) install.packages("ape")
if (!require("TreeDist")) install.packages("TreeDist")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("Quartet")) install.packages("Quartet")
if (!require("phangorn")) install.packages("phangorn")
if (!require("tidyr")) install.packages("tidyr")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

message("Loading in the tree files...\n")
# true tree (cellcoal) and inferred_trees (beast_tree and cellphy_tree)
true_tree <- read.tree(TRUETREE_PATH)
beast_tree <- read.nexus(BEASTtree_PATH)
cellphy_tree <- read.tree(CELLPHYtree_PATH)

message("Filtering out any uncommon tip labels from trees...\n")
# Removing a "cell" prefix from CellCoal tip lables 
# Removing a "tumecell" prefix from CellPhy/ BEAST tip labels

true_tree$tip.label <- gsub("cell", "", true_tree$tip.label)
beast_tree$tip.label <- gsub("tumcell", "",beast_tree$tip.label)
cellphy_tree$tip.label <- gsub("tumcell", "", cellphy_tree$tip.label)

#Finding common tips between trees and ensuring they only include those common tips
common_tips <- Reduce(intersect, list(
  true_tree$tip.label, 
  cellphy_tree$tip.label, 
  beast_tree$tip.label
))

true_tree <- keep.tip(true_tree, common_tips)
beast_tree <- keep.tip(beast_tree, common_tips)
cellphy_tree <- keep.tip(cellphy_tree, common_tips)

message("Calculating Topological Metrics between true and inferred trees...\n")
# Calculating different Topological Comparisons between true tree and inferred trees (How accurate are the splits?)

# Normalised RF Distances, highly senstive and can be seen to waiver dramatically 

RF_beast <- TreeDist::RobinsonFoulds(true_tree, beast_tree, normalize = TRUE)
RF_cellphy <- TreeDist::RobinsonFoulds(unroot(true_tree), cellphy_tree, normalize = TRUE)

# Generalized RF (Clustering Information Distance is more sensitive)
# GRF is more robust and less susceptible to the "lack of resolution" or small changes causing a massive jump in score

CID_beast   <- ClusteringInfoDistance(beast_tree, true_tree, normalize = TRUE)
CID_cellphy <- ClusteringInfoDistance(cellphy_tree, unroot(true_tree), normalize = TRUE)

# Quartet Divergence (Best for single-cell local accuracy)
# QuartetStatus, calculates how many quartets are identical, differ, or are unresolved
# QuartetDivergence, takes the status and calculates the symmetric quartet divergence (0 = identical, 1 = maximum possible difference)

status_beast <- QuartetStatus(beast_tree, true_tree)
status_cellphy <- QuartetStatus(cellphy_tree, unroot(true_tree))

QD_beast   <- QuartetDivergence(status_beast, similarity = FALSE)
QD_cellphy <- QuartetDivergence(status_cellphy, similarity = FALSE)

message("Generating and storing Metrics data frame in a .csv file within similarity_tests folder\n")
# Creating Data Frame of Metrics and saving to .csv for future reference 
results <- data.frame(
  Model = c("BEAST (Bayesian)", "CellPhy (ML)"),
  RF = c(RF_beast, RF_cellphy),
  CID = c(CID_beast, CID_cellphy),
  QD = c(QD_beast, QD_cellphy)
)

save_path <- file.path(OUTPUT_PATH, "ComparisonMetric_Results.csv")
write.csv(results, file = save_path, row.names = FALSE)

message("Generating Clustered Barplot...\n")
#Creating Clustered BarPlot 

# Pivoting to long format to seperate columns based on Model used, Combine RF,CID and QD metrics onto the x-axis and Normalised values on y axis 
results_long <- results %>%
  pivot_longer(cols = -Model, 
               names_to = "Metric", 
               values_to = "Value")

# Generating and formatting clustered barplot of comparison metric results 
ClusBarplot <- ggplot(results_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", width=0.7) +  #"dodge" ensures bars are beside each other and 0.7 width of the bars 
    scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +  # Adds padding to the left and right of plot
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +   
  theme(
    text = element_text(family = "serif"), # Setting the font family
    plot.title = element_text(size = 20, hjust = 0.7, face = "bold", color = "black"), #size of title, hjust: position of title across page
    axis.line = element_line(color = "black", linewidth = 1.5), # Makes axis lines black and thicker
    axis.ticks = element_line(color = "black", linewidth = 1),  # Makes axis ticks black and thicker
    axis.title.x = element_text(size = 14, face = "bold", color = "black"), #Formating the axis titles
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12), #Formatting the text for Metric names
    axis.text.y = element_text(size = 12), #Formatting the text for Values
    legend.title = element_text(size = 16, face = "bold"), # Formatting Legend title and text
    legend.text = element_text(size = 14),  
    plot.margin = unit(c(2, 1.5, 1, 1.5), "cm"), # Adjusting Top, Right, Bottom, Left margins 
    plot.background = element_rect(fill = "gray90"), #Setting plot background color 
    panel.background = element_rect(fill = "gray90"),
  ) +
  labs(title = "Comparing Phylogenetic Fidelity of Cellphy and BEAST2 in
   Reconstructing CellCoal-Simulated Cancer Evolution Trees\n",
       y = "Normalised Distance Value\n",
       x = "\nMetric Type") +
       scale_fill_brewer(palette = "Set2") 

message("Saving Clustered Barplot in .pdf format within similarity_tests folder\n")
save_path <- file.path(OUTPUT_PATH, "MetricComparison_plot.pdf")
ggsave(filename = save_path, plot = ClusBarplot, width = 12, height = 10, units = "in")

message("Done!")


----
