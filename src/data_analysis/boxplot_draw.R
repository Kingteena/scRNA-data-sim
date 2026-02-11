## Final Step of Simulation Workflow
## R script summarises the accuracy of phylogentic tree inferences from CellPhy and BEAST through key comparison  metrics
## Boxplot (.pdf) and Table (.csv) Stats Summary (.txt) are generated for each replicates set of three trees inputted from script caller


# trailingOnly = TRUE: ignores all the background "system" info and only focuses on specific arguments you typed in Bash in the caller 
# args:  Creates a vector of each argument as a numbered item 

args <- commandArgs(trailingOnly = TRUE)
# Rscript $BOXPLOT_BEAST_SCRIPT "$true_tree_path" "$beast_trees_path" "cellphy_trees_path" "$output_path" "$RUN_NUM"

True_trees_path    <- args[1]
Beast_trees_path   <- args[2]
Cellphy_trees_path <- args[3]
output_path        <- args[4]
run_number         <- args[5]
 
# Load Libraries
suppressPackageStartupMessages({
  library(ape)
  library(TreeDist)
  library(ggplot2)
  library(gtools)
  library(Quartet)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
})

# Custom Functions for Tree Processing 

# process_tree : Function to process the inputted trees; matching tip labels and normalising tree height 
process_tree <- function(tree) {
  # Strip everything except digits to match 'tumcell_0001' and 'cell0001'
  new_labels <- gsub("[^0-9]", "", tree$tip.label)

  # Remove invisible whitespace/carriage returns
  new_labels <- trimws(new_labels)

  # Rename the outgroup (empty strings) to "OUTGROUP"
  new_labels[new_labels == ""] <- "OUTGROUP"
  tree$tip.label <- new_labels

 # Normalise branch lengths to total height of 1.0
  max_height <- max(ape::node.depth.edgelength(tree), na.rm = TRUE)
  if(max_height > 0) tree$edge.length <- tree$edge.length / max_height  

  return(tree) 
}

# Metrics Calculation Function
# 5 key metrics (RF,JRF,QD,CID,Path_Disr) for assessing accuracy of the phylogentic inference

calculate_metrics <- function(clean_true, clean_inferred, method_name) {

  # Syncing tip labels to avoid "length mismatch" errors 
  # Saving into a "current_"tree, to avoid corrupting the original tree files
  # 'intersect' automatically ignores the 21st Outgroup cell because it exists in CellCoal but NOT in the BEAST CCD0  or Cellphy trees
  common <- intersect(clean_true$tip.label, clean_inferred$tip.label)
  current_true <- ape::keep.tip(clean_true, common)
  current_inferred <- ape::keep.tip(clean_inferred, common)

  # Unrooting both trees for topological metrics; RF/CID/QD/JRF are best calculated on unrooted trees
  current_true <- ape::unroot(current_true)
  current_inferred <- ape::unroot(current_inferred)

  # Calculating Key Metrics 
  RF_value  <- TreeDist::RobinsonFoulds(current_true, current_inferred, normalize = TRUE)

  JRF_value <- TreeDist::JaccardRobinsonFoulds(current_true, current_inferred, k = 1.1, normalize = TRUE)

  CID_value <- TreeDist::ClusteringInfoDistance(current_true, current_inferred, normalize = TRUE)

  status_value <- Quartet::QuartetStatus(current_true, current_inferred)
  QD_value     <- Quartet::QuartetDivergence(status_value, similarity = FALSE)
  
  base_path <- TreeDist::PathDist(current_true, current_inferred)

  # Normalization for Path Distance : Error is relative to the size of the tree 
  n_tips <- length(common)  # total number of taxa (tips) in the tree
  max_path_possible <- sqrt(choose(n_tips, 2) * (n_tips - 1)^2)  #Root mean squared error of path lengths (Steel and Penny Distance)
  Path_dist_norm <- base_path / max_path_possible  # Standardisation 

  return(data.frame(
    Method = method_name,
    RF = RF_value, 
    JRF = JRF_value, 
    Quartet = QD_value, 
    CID = CID_value, 
    PathDist = Path_dist_norm
  ))
}

# Main Script
# 1. Loading in Tree Files 

message("Loading Tree files for BEAST and CellPhy comparison...")

# Enlisitng "gtools::mixedsort" to ensure tree files match (e.g 0001 matches 0001), when comparing Beast and Cellcoal trees
# for true_trees, trees\\ = matching with the text 'trees', [0-9]+ = matches one or more digits, $ = the file has to end there, preventing log files

truetree_files    <- gtools::mixedsort(list.files(True_trees_path, pattern = "trees\\.[0-9]+$", full.names = TRUE))
beasttree_files   <- gtools::mixedsort(list.files(Beast_trees_path, pattern = ".*\\.tree$", full.names = TRUE))
cellphytree_files <- gtools::mixedsort(list.files(Cellphy_trees_path, pattern = ".*\\.vcf\\.raxml\\.bestTree$", recursive = TRUE,  full.names = TRUE))

#Starting an empty list for R to place data after every run 
results_list <- list()


message("Looping through comparing trees found from loading in directories\n")

# Enlisting "seq_along" instead of 1:length(x), as if the folder is 0, it doesn't start the loop 
for (i in seq_along(truetree_files)) {

  Tree_ID <- gsub(".*\\.([0-9]+)$", "\\1", basename(truetree_files[i]))

  #Reading in Ground Truth Tree (CellCoal)
  # With some low mutation rate (e.g run_1), some basepairs don't acquire enough mutations to define a tree. 
  # If CellCoal doesn't find any mutations, it instead produces an "empty" tree or a tree that ape cannot parse, thus if statement is added to check for this 
 
    true_tree <- tryCatch({
    tree <- ape::read.tree(truetree_files[i])
    if(is.null(tree) || class(tree) != "phylo") stop("Empty")
    process_tree(tree)
  }, error = function(e) return(NULL))

  if (is.null(true_tree)) {
    message("!! Skipping ID ", Tree_ID, ": True tree file is empty or unreadable.")
    next
  }

  # Reading in Inference Tree (BEAST)
  # Within in the loop; reading in trees, processing, calculating metrics and adding row to results_list
  # TryCatch: is a safety net which tries Nexus read in first, then Newick incase of error
  # Added safety measure in case of "empty" trees from Cellcoal 
  # if(is.null(beast_raw) || class(beast_raw) != "phylo") stop("Empty") = if tree file has an error message, it triggers the stop command, prints message and assigns NULL to beast_tree

  if (i <= length(beasttree_files)) {
    beast_tree <- tryCatch({
      beast_raw <- tryCatch({ape::read.nexus(beasttree_files[i])}, 
                         error = function(e) ape::read.tree(beasttree_files[i]))
      if(is.null(beast_raw) || class(beast_raw) != "phylo") stop("Empty")
      process_tree(beast_raw)
    }, error = function(e) return(NULL))

    if (!is.null(beast_tree)) {
    beast_metrics <- calculate_metrics(true_tree, beast_tree, "BEAST")
    beast_metrics$Tree_ID <- Tree_ID
    results_list[[length(results_list) + 1]] <- beast_metrics
    } else {
      message (" Warning: BEAST tree ", Tree_ID, "is missing. ")
    }
  }

  # Reading in Inference Tree (CellPhy)
  # Within in the loop; reading in trees, processing, calculating metrics and adding row to results_list
  # Added safety measure in case of "empty" trees from Cellcoal 

 if (i <= length(cellphytree_files)) {
    cellphy_tree <- tryCatch({
      cellphy_raw <- ape::read.tree(cellphytree_files[i])
      if(is.null(cellphy_raw) || class(cellphy_raw) != "phylo") stop("Empty")
      process_tree(cellphy_raw)
    }, error = function(e) return(NULL))

  if (!is.null(cellphy_tree)){
      cellphy_metrics <- calculate_metrics(true_tree, cellphy_tree, "CellPhy")
      cellphy_metrics$Tree_ID <- Tree_ID
      results_list[[length(results_list) + 1]] <- cellphy_metrics
  } else{
    message("Warning: Cellphy Tree", Tree_ID, "is missing.")
  }
}
  # The modulo operator (%%), divides iteration by 20, if the remainder is 0, then it prints the message (can be adjust for preference)
  # Keeps terminal easy to read 
  if(i %% 20 == 0) message("Completed analysis for Tree pair: ", Tree_ID)

}


# Output Data and Statistics Calculations

#Checking if there were valid tree pairs found within the run before combining 
if (length(results_list) == 0) {
  stop("STOP: No valid tree pairs were found in the run, Please check simulated logs in output folder! ")
}
# Combining rows from results list together to form metrics_combined file 
metrics_combined <- dplyr::bind_rows(results_list)

message("Performing statistical comparisons...")

# Calculating Averages and Standard Deviation Values for Cellphy and BEAST 
stats_summary <- metrics_combined %>%
  tidyr::pivot_longer(cols = c(RF, JRF, Quartet, CID, PathDist), names_to = "Metric", values_to = "Value") %>%
  group_by(Metric, Method) %>%
    dplyr::summarise(
    Average = mean(Value, na.rm = TRUE),
    StdDev  = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = Method, values_from = c(Average, StdDev))


# We pivot wide so that BEAST and CellPhy values for the same Tree_ID are on one row
# Employing AND Performing Wilcoxon Signed Rank Paired Test
# In high SNR runs, many replicates will produce the exact distance score for Cellphy and Beast, leading to wilcoxon test not being able to "rank" a zero difference; indicating performance convengence 

stats_wilcoxon <- metrics_combined %>%
  tidyr::pivot_longer(cols = c(RF, JRF, Quartet, CID, PathDist), names_to = "Metric", values_to = "Val") %>%
  tidyr::pivot_wider(names_from = Method, values_from = Val) %>%
  group_by(Metric) %>%
  dplyr::summarise(
     p_value = tryCatch({
      wilcox.test(BEAST, CellPhy, paired = TRUE, exact = FALSE)$p.value  #exact = FALSE, prevents terminal cluttering
    }, error = function(e) return(NA)), 
    .groups = "drop"
  )

# Combining stats and p.values into a final overall results table 
# left_join: using 'Metric' column as the key to merge Mean/SD table with the p-value table

final_stats <- dplyr::left_join(stats_summary, stats_wilcoxon, by = "Metric")

# Sink function drives output from R to the external file 'connection'
# Writing to the Final Summary File 

sink(file.path(output_path, "statistical_summary.txt"))
cat("Phylogenetic Accuracy Summary Report: ", run_number, "\n")
cat("=======================================================================================================\n")
cat("Metrics represent Normalized Distance to Cellcoal True Tree (Aiming to be as close to Zero as possible)\n")
cat("Statistical Test employed is the Paired Wilcoxon Signed-Rank (BEAST vs CellPhy)\n")
cat("======================================================================================================\n\n")

# Limiting decimal points for clarity and ease of read
print(as.data.frame(final_stats), digits = 4)

cat("\nNote: p_value = NA indicates identical results across all replicates.\n")
sink()

# Preparing Dataframe for Boxplot Plotting

# Combining into a single 'Long' dataframe for ggplot2 and statistics calculations 
metrics_long <- metrics_combined %>%
  tidyr::pivot_longer(cols = c(RF, JRF, Quartet, CID, PathDist), 
                     names_to = "Metric", values_to = "Distance") 

# Saving combined dataframe into .csv file format and storing in output path
write.csv(metrics_long, file.path(output_path, "metric_calculations_combined.csv"), row.names = FALSE)

#Ensuring that the metrics on the x-axis of the boxplots, are in a particular order (change per preference)
metrics_long$Metric <- factor(metrics_long$Metric, levels=c("RF", "JRF", "CID", "Quartet", "PathDist"))

#Calculating the upper_limit on the y-axis for this plot, in order to adjust the plot accordingly 
upper_limit <- max(metrics_long$Distance, na.rm = TRUE) * 1.1

# Generating the Boxplot for Calculated Metrics 

metric_boxplot <- ggplot(metrics_long, aes(x = Metric, y = Distance, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.9)) +

  # Geom_jitter - spread the 0.0 values so we can see the density
  geom_jitter(aes(group = Method), position = position_dodge(width = 0.9), 
              alpha = 0.3, size = 1.2, color = "midnightblue", na.rm = TRUE) +

  # geom_hline - Acts as a Reference line at 0, the goal for all metrics to be 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  
  coord_cartesian(ylim = c(0, upper_limit))+
  theme_minimal() +
  theme(
    text = element_text(family = "serif"), # Setting the font family
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold", color = "black"), #size of title, hjust: position of title across page
    axis.title.x = element_text(size = 14, face = "bold", color = "black"), #Formating the axis titles
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12), #Formatting the text for Metric names
    axis.text.y = element_text(size = 12), #Formatting the text for Values
    legend.title = element_text(size = 14, face = "bold"), # Formatting Legend title and text
    legend.text = element_text(size = 14),  
    plot.margin = unit(c(2, 2, 1, 2), "cm"), # Adjusting Top, Right, Bottom, Left margins 
  ) + 
  #scale_fill_manual(values = c("BEAST" = "#66c2a5", "CellPhy" = "#fc8d62")) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Comparing Phylogenetic Accuracy of BEAST (CCD0) and CellPhy (ML) in 
    Reconstructing CellCoal-Simulated Cancer Evolution Trees\n",
       subtitle = "Data pooled from 100 Replicates ",
       y = "Normalized Distance\n", # value of 0 = identical topologies 
       x = "\nAccuracy Metric Type") +
      theme(legend.position = "top")
       


ggsave(file.path(output_path, "Metric_Comparison_Plot.pdf"), metric_boxplot, width = 11, height = 7)

message("Analysis complete. Results saved to: ", output_path)

