#!/usr/bin/env Rscript

# 05-nonmalignant_clustering_and_qc_plots.R
library(Seurat)

non_malignant <- readRDS("data/processed/non_malignant_processed.rds")

# Determine PCs, neighbors and clusters (parameters tunable)
non_malignant <- FindNeighbors(non_malignant, dims = 1:15)
non_malignant <- FindClusters(non_malignant, resolution = 0.1, random.seed = 40)

# Save cluster annotations
saveRDS(non_malignant, file = "data/processed/non_malignant_clustered.rds")

# Create and save some diagnostic plots programmatically (png files)
plots_dir <- "figures"
dir.create(plots_dir, showWarnings = FALSE)

png(file.path(plots_dir, "nonmalignant_umap_celltype.png"), width = 800, height = 800)
print(DimPlot(non_malignant, reduction = "umap", group.by = "cell_type", pt.size = 1.2))
dev.off()

png(file.path(plots_dir, "nonmalignant_umap_tumorid.png"), width = 800, height = 800)
print(DimPlot(non_malignant, reduction = "umap", group.by = "tumor_id", pt.size = 1.2))
dev.off()

message("Saved clustering results and UMAP plots in figures/")
