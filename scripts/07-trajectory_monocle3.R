#!/usr/bin/env Rscript

# 07-trajectory_monocle3.R
# This script converts a Seurat object (malignant subset) to a Monocle3 CellDataSet and runs
# a basic trajectory. Requires SeuratWrappers and monocle3 installed.

library(Seurat)
library(SeuratWrappers)
library(monocle3)

malignant <- readRDS("data/processed/malignant_processed.rds")

# Example: choose Mel79 (as in your Rmd)
if (!"tumor_id" %in% colnames(malignant@meta.data)) stop("tumor_id not found in metadata")
mel79 <- subset(malignant, subset = tumor_id == "Mel79")

cds <- as.cell_data_set(mel79)

# ensure gene_short_name is present
gene_names <- rownames(mel79[[]])
if (is.null(cds@rowRanges@elementMetadata$gene_short_name)) {
  cds@rowRanges@elementMetadata$gene_short_name <- gene_names
}

cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 50, method = 'LSI')
cds <- reduce_dimension(cds, reduction_method = 'UMAP', preprocess_method = 'LSI')
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# pick a root programmatically
# Example: choose a cell with minimal pseudotime after learn_graph
# cds <- order_cells(cds)

saveRDS(cds, file = "data/processed/mel79_monocle_cds.rds")
message("Saved Monocle3 CDS to data/processed/mel79_monocle_cds.rds")
