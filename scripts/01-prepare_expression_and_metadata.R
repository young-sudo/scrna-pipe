#!/usr/bin/env Rscript

# 01-prepare_expression_and_metadata.R
library(dplyr)
library(Seurat)

infile <- "data/processed/raw_tirosh_table.rds"
raw <- readRDS(infile)

# The input file appears to have 3 header rows and then gene x cell matrix
# Adjust according to your file structure (this follows your Rmd)
meta_rows <- raw[1:3, , drop = FALSE]
gene_table <- raw[-(1:3), , drop = FALSE]

# Make unique gene names from first column
gene_names <- make.unique(gene_table[[1]], sep = ".")
rownames(gene_table) <- gene_names
expr <- gene_table[, -1]

# Ensure numeric matrix
expr_mat <- as.matrix(sapply(expr, as.numeric))
rownames(expr_mat) <- gene_names

# Transpose and build metadata
meta <- as.data.frame(t(meta_rows))
rownames(meta) <- meta[[1]]
meta <- meta[, -1, drop = FALSE]

# Basic cleaning of suspicious columns (keep your original checks)
cell_names <- colnames(expr_mat)
extracted_info <- sapply(cell_names, function(x) paste0("Mel", substr(x, 3, 4)))
incorrect_entries <- which(extracted_info %in% c("Mel2_", "Melni"))
if (length(incorrect_entries) > 0) {
  expr_mat <- expr_mat[, -incorrect_entries]
  meta <- meta[setdiff(rownames(meta), names(incorrect_entries)), , drop = FALSE]
}

# Save expression and metadata
saveRDS(list(counts = expr_mat, meta = meta), file = "data/processed/expr_and_meta.rds")
message("Saved expression and metadata to data/processed/expr_and_meta.rds")
