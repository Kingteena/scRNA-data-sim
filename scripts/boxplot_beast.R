## Final Step of Simulation Workflow

args <- commandArgs(trailingOnly = TRUE)
# Rscript $BOXPLOT_BEAST_SCRIPT "$true_tree_path" "$beast_trees_path" "$output_path" "$RUN_NUM"

True_tree_path <- args[1]
Beast_trees_path <- args[2]
output_path <- args[3]
run_number <- args[4]

#True_tree_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/cellcoal/run_3/trees_dir/0001.tree"
#Beast_trees_path="/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/beast/run_3/tree/1.tree"
#output_path = "/hpcfs/groups/phoenix-hpc-gavryushkina/simulation/output_sahana/beast/run_3/metric_calcs"

message("Loading in any required libraries into current workspace")
# Loading in any required libraries into current workspace 

if (!require("ape")) install.packages("ape")
if (!require("TreeDist")) install.packages("TreeDist")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("Quartet")) install.packages("Quartet")
if (!require("phangorn")) install.packages("phangorn")
if (!require("tidyr")) install.packages("tidyr")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")


# Function to process the inputted the trees; matching tip labels, normalising tree height and unrooting tree 
process_tree <- function(tree, type = "true") {

    message("Filtering out tip label suffix's from trees...\n")
  # Removing tumecell (BEAST) or cell (CellCoal) prefix from Tip labels 
  
     tree$tip.label <- gsub("tumcell", "", tree$tip.label, ignore.case = TRUE) 
     tree$tip.label <- gsub("cell", "", tree$tip.label, ignore.case = TRUE) 

    message("Normalising Total tree depth...\n")
  # Normalizing Tree Heights by scales total tree depth to maximum 1
  # KF (Branch Score) can then be comparable across different simulations
  # Doesn't affect QD, RF and CID metrics as they are Toplogy based similairty 
  max_height <- max(ape::node.depth.edgelength(tree))
   if(max_height > 0) tree$edge.length <- tree$edge.length / max_height
  
  message("Returning unrooted tree...\n")
  # Unrooting trees for topologival and distance metrics
  # RF/CID/QD/KF are best calculated on unrooted trees
  return(unroot(tree))
}

# Metrics Calculation Function
# 4 key metrics (QD,RF,KF,CID) for assessing accuracy of phylogentic inference
calculate_metrics <- function(clean_true, clean_beast) {

 # Ensuring Original tree remains untouched for the next iteration
  current_true  <- clean_true
  current_beast <- clean_beast

 # Syncing tips to avoid "length mismatch" errors 
  common <- intersect(current_true$tip.label, current_beast$tip.label)
  current_true <- keep.tip(current_true, common)
  current_beast <- keep.tip(current_beast, common)


  message("Calculating Similarity Metrics...\n")
  # Normalised RF Distances, highly senstive and can be seen to waiver dramatically 
  # Using package from TreeDist due to normalisation ability

  RF_beast <- TreeDist::RobinsonFoulds(current_true, current_beast, normalize = TRUE)

  # Quartet Divergence (Best for single-cell local accuracy)
  # QuartetStatus, calculates how many quartets are identical, differ, or are unresolved
  # QuartetDivergence, takes the status and calculates the symmetric quartet divergence (0 = identical, 1 = maximum possible difference)

  status_beast <- QuartetStatus(current_beast, current_true)
  QD_beast <- QuartetDivergence(status_beast, similarity = FALSE)
  
  # Clustering Information Distance (CID); How much "information" the two trees share regarding their grouped taxa 
  # CID is less susceptible to the "lack of resolution" or small changes in leafs causing a massive jump in score

  CID_beast <- TreeDist::ClusteringInfoDistance(current_beast, current_true, normalize = TRUE)

  # D. Kuhner-Felsenstein (KF) Value (Branch Score Distance)
  # Similarity of trees in terms of the evolutionary distance; by comparing them by their topology and branch lengths (evolutionary time/divergence)

  KF_beast <- phangorn::KF.dist(current_beast, current_true)
  #total_len <- sum(current_true$edge.length) + sum(current_beast$edge.length)
  #KF_beast <- if(total_len > 0) raw_kf_beast/total_len else NA  #avoiding dividing by zero 

  message("Returning Data Frame of calculated similarity metrics..\n")

  return(data.frame(RF = RF_beast, Quartet = QD_beast, CID = CID_beast, KF = KF_beast))
}

message("Loading in and Processing CellCoal origin tree file...\n")

true_tree <- read.tree(True_tree_path)
clean_true <- process_tree(true_tree)

message("Loading in the BEAST tree files...\n")

beasttree_files <- list.files(Beast_trees_path,
                   pattern = ".*\\.tree$",  
                   full.names = TRUE)

#Starting an empty list for R to place data after every run 
results_list <- list()

message("Looping through trees found from loading directories\n")

for (i in 1:length(beasttree_files)) {
  # do we too handle potential multiPhylo from BEAST?
   #raw_beast <- read.nexus(beasttree_files[i])
  #if(inherits(raw_beast, "multiPhylo")) raw_beast <- raw_beast[[length(raw_beast)]]
  
  beast_tree <- read.nexus(beasttree_files[i])
  clean_beast <- process_tree(beast_tree)  # processing each BEAST file per iteration 
  
  # Verifying tip labels match
  # Otherwise metrics default to 1.0 (max), indicating error in process_tree function
  # setequalcheck is two objects contain exactly the same unique elements

  if(!setequal(clean_true$tip.label, clean_beast$tip.label)) {
      message("Note: Tree ", i, " has different tip counts. Pruning to common subset...")
  }

  similarity_metrics <- calculate_metrics(clean_true, clean_beast)  #running function for each Beast tree with Origin Tree 
  similarity_metrics$Tree_Number <- paste0("Tree", i) #Adding new column "Trees" to results, then "paste0" glues "Tree" to the current interation
  message("Completed metric calculations for Tree ", i, ": ", beasttree_files[i]) # interation number : tree file name
  results_list[[i]] <- similarity_metrics  # Prevents overwriting 
  }

# Combining into a single 'Long' dataframe for ggplot2

metrics_combined <- dplyr::bind_rows(results_list) %>%
pivot_longer(cols = c(RF, Quartet, CID, KF), names_to = "Metric", values_to = "Distance")


message("Saving Calculated metrics into a .csv file\n")

# Saving combined dataframe into .csv file and storing in output path

save_path <- file.path(output_path, "Metric_Calculations.csv")
write.csv(metrics_combined, file = save_path, row.names = FALSE)

# Finding Required Adjustment height for boxplot output 
# Finding the highest value in your specific column
max_val <- max(metrics_combined$Distance)  #na.rm = TRUE

# Adding a 10% buffer 
upper_limit <- max_val * 1.1

# Generating Boxplot

message("Now Generating Boxplot of Calculated metrics\n")

metric_boxplot<- ggplot(metrics_combined, aes(x = Metric, y = Distance, fill = Metric)) +
   geom_boxplot(outlier.shape = NA) +
  #coord_cartesian(ylim = c(0, upper_limit)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +
  ylim(0, 1) +
  theme_minimal() +
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
       subtitle = "Data pooled from 3 replicates",
       y = "Normalized Distance\n", # value of 0 = identical topologies 
       x = "\nMetric Type") +
      scale_fill_brewer(palette = "Set2") 
       
#Saving boxplot output, into output folder 

message("Saving Boxplot in .pdf format\n")
save_path <- file.path(output_path, "Metric_Calculations_plot.pdf")
ggsave(filename = save_path, plot = metric_boxplot, width = 12, height = 20, units = "in")

sprintf("The .csv data frame and .pdf boxplot are found in run_%s/beast/metric_calcs folder", run_number)
message("Done!")


  


