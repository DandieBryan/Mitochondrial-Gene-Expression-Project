### MADE BY BRYAN LE, JULY 2024

# This part of the pipeline compares and visualizes the gene expression of mitochondrial 
# genes calculated by two different programs, 'RNA-SeQC' and 'Salmon'



##########################################################################################
# Load libraries
.libPaths("INSERT PATH TO R LIBRARY")
library(ggplot2)
library(reshape2)
library(gridExtra)  
library(grid)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
rnaseqc_ge_file <- args[1]
salmon_ge_file <- args[2]
output_directory <- args[3]
reference_genes_file <- args[4]

# Read the expression data
data1 <- read.table(rnaseqc_ge_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
data2 <- read.table(salmon_ge_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Reshape data to long format
data1_long <- melt(data1, id.vars = "GeneID", variable.name = "SampleID", value.name = "TPM_x")
data2_long <- melt(data2, id.vars = "GeneID", variable.name = "SampleID", value.name = "TPM_y")

# Read the gene of interest file
genes_of_interest <- read.csv(reference_genes_file, stringsAsFactors = FALSE)

# Filter for genes with transcript type "protein_coding"
filtered_genes <- subset(genes_of_interest, transcript_type == "Mt_tRNA")

# Extract the list of gene_ids for the genes of interest
genes_to_include <- filtered_genes$gene_id

# Filter data to include only the genes of interest
data1_mt_genes <- subset(data1_long, GeneID %in% genes_to_include)
data2_mt_genes <- subset(data2_long, GeneID %in% genes_to_include)

# Find common genes
common_genes <- intersect(unique(data1_mt_genes$GeneID), unique(data2_mt_genes$GeneID))

# Filter data to include only common genes
data1_filtered <- subset(data1_mt_genes, GeneID %in% common_genes)
data2_filtered <- subset(data2_mt_genes, GeneID %in% common_genes)

# Find common samples
common_samples <- intersect(data1_filtered$SampleID, data2_filtered$SampleID)

# Further filter data to include only common samples
data1_final <- subset(data1_filtered, SampleID %in% common_samples)
data2_final <- subset(data2_filtered, SampleID %in% common_samples)

# Merge the filtered data on SampleID and GeneID
merged_data <- merge(data1_final, data2_final, by = c("SampleID", "GeneID"))

# Apply log transformation with a small offset (if necessary)
merged_data$TPM_x <- log10(merged_data$TPM_x + 1)
merged_data$TPM_y <- log10(merged_data$TPM_y + 1)

# List to store plots
plots <- list()

# Generate scatter plots for each gene
for (gene_id in unique(merged_data$GeneID)) {
  # Get the gene_name corresponding to the gene_id
  gene_info <- filtered_genes[filtered_genes$gene_id == gene_id, ]
  gene_name <- gene_info$gene_name[1]  # Extract the gene_name (assuming there's only one match)
  
  # Filter data for the current gene
  gene_data <- subset(merged_data, GeneID == gene_id)
  
  # Perform linear regression on the log-transformed data 
  model <- lm(TPM_y ~ TPM_x, data = gene_data)
  summary(model)
  
  # Create the scatterplot for the current gene
  p <- ggplot(gene_data, aes(x = TPM_x, y = TPM_y)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste("Scatterplot of", gene_name, "Expression"),
         x = "Log10(Gene Expression in TPM from RNA-SeQC + 1)",
         y = "Log10(Gene Expression in TPM from Salmon + 1)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),  # Increase title font size
      axis.title = element_text(size = 18),  # Increase axis titles font size
      axis.text = element_text(size = 14),  # Increase axis text font size
      plot.margin = unit(c(1, 1, 1, 1), "cm")  # Add margins to each plot
      )
  
  # Store the plot in the list
  plots[[length(plots) + 1]] <- p
}

# Create a 2x7 grid of plots and save it
n <- length(plots)
ncol <- 3
nrow <- ceiling(n / ncol)

# Define the title using grid::textGrob with additional margin at the top
title <- grid::textGrob("Expression of Mitochondrial tRNA Genes Estimated by RNA-seQC versus Salmon",
                        gp = gpar(fontsize = 48, fontface = "bold"),
                        y = unit(1, "npc") - unit(3.0, "cm"))

# Arrange and save the plots
png(file.path(output_directory, "MT_tRNA_GE_scatterplots.png"), width = 16 * ncol, height = 9 * nrow, units = "in", res = 300)
grid.arrange(grobs = plots, ncol = ncol, nrow = nrow, top = title, padding = unit(8, "lines"))  # Add title and adjust padding
dev.off()
