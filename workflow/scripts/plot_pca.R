#!/usr/bin/env Rscript
#
# Generate PCA plot from count matrix
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
  library(ggplot2)
})

# Read count matrix
cat("Reading count matrix...\n")
counts <- read.table(snakemake@input$counts, header = TRUE, row.names = 1, sep = "\t")
cat("Count matrix dimensions:", dim(counts), "\n")

# Read sample metadata
cat("Reading sample metadata...\n")
metadata <- read.table(snakemake@input$metadata, header = TRUE, sep = "\t")
cat("Sample metadata:\n")
print(metadata)

# Filter low-count genes
cat("\nFiltering low-count genes...\n")
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# Select top variable genes for PCA
cat("\nSelecting top variable genes...\n")
n_genes <- min(snakemake@params$n_top_genes, nrow(counts_filtered))
gene_vars <- apply(counts_filtered, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_genes])
counts_top <- counts_filtered[top_genes, ]
cat("Using top", n_genes, "most variable genes\n")

# Log-transform counts
counts_log <- log2(counts_top + 1)

# Perform PCA
cat("\nPerforming PCA...\n")
pca_result <- prcomp(t(counts_log), scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create data frame for plotting
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$condition <- metadata$condition[match(pca_df$sample, metadata$sample)]

cat("\nPCA variance explained:\n")
cat("PC1:", round(var_explained[1], 2), "%\n")
cat("PC2:", round(var_explained[2], 2), "%\n")

# Create PCA plot
cat("\nGenerating PCA plot...\n")
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1, size = 3) +
  labs(
    title = "PCA of RNA-Seq Samples",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Condition"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Save plot
ggsave(snakemake@output[[1]], plot = p, width = 8, height = 6)
cat("PCA plot saved to", snakemake@output[[1]], "\n")

# Close log connections
sink(type = "message")
sink(type = "output")
close(log_file)
