#!/usr/bin/env Rscript
#
# Generate sample correlation heatmap
#
# This script is called by Snakemake and has access to:
# - snakemake@input: list of input files
# - snakemake@output: list of output files
# - snakemake@params: parameters from the rule
# - snakemake@log: log file path

# Redirect output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# Load required libraries
suppressPackageStartupMessages({
  library(pheatmap)
})

# Read count matrix
cat("Reading count matrix...\n")
counts <- read.table(snakemake@input$counts, header = TRUE, row.names = 1, sep = "\t")
cat("Count matrix dimensions:", dim(counts), "\n")

# Read sample metadata
cat("Reading sample metadata...\n")
metadata <- read.table(snakemake@input$metadata, header = TRUE, sep = "\t", row.names = 1)

# Filter low-count genes
cat("\nFiltering low-count genes...\n")
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# Log-transform counts
counts_log <- log2(counts_filtered + 1)

# Calculate correlation matrix
cat("\nCalculating sample correlations...\n")
cor_method <- snakemake@params$method
cor_matrix <- cor(counts_log, method = cor_method)

cat("Correlation matrix:\n")
print(round(cor_matrix, 3))

# Prepare annotation
annotation_col <- data.frame(
  Condition = metadata$condition,
  row.names = rownames(metadata)
)

# Generate heatmap
cat("\nGenerating correlation heatmap...\n")
pdf(snakemake@output[[1]], width = 8, height = 7)

pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(0.8, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.3f",
  fontsize = 10,
  main = paste0("Sample Correlation Heatmap (", cor_method, ")"),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)

dev.off()

cat("Correlation heatmap saved to", snakemake@output[[1]], "\n")

# Close log connections
sink(type = "message")
sink(type = "output")
close(log_file)
