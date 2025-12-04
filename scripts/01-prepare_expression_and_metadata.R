#!/usr/bin/env Rscript

# 01-prepare_expression_and_metadata.R
library(dplyr)
library(Seurat)
library(optparse)

infile <- "results/processed/raw_table.rds"
out_dir <- "results/processed/expressions"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=infile, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
infile <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- readRDS(infile)

# The input file appears to have 3 header rows and then gene x cell matrix
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

outfile <- file.path(out_dir, "expr_and_meta.rds")
# Save expression and metadata
saveRDS(list(counts = expr_mat, meta = meta), file = outfile)
message("Saved expression and metadata to ", outfile)
