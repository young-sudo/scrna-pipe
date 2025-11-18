#!/usr/bin/env Rscript

# 03-filter_and_split_malignant.R
library(Seurat)

seurat_obj <- readRDS("data/processed/seurat_before_filtering.rds")

# If you want to apply threshold-based filtering, add here (example thresholds)
# seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# The dataset contains a metadata column used in your Rmd:
# malignant.1.no.2.yes.0.unresolved.
if (!"malignant.1.no.2.yes.0.unresolved." %in% colnames(seurat_obj@meta.data)) {
  stop("Expected metadata column malignant.1.no.2.yes.0.unresolved. not found")
}

non_malignant <- subset(seurat_obj, subset = malignant.1.no.2.yes.0.unresolved. == 1)
malignant <- subset(seurat_obj, subset = malignant.1.no.2.yes.0.unresolved. == 2)

saveRDS(list(non_malignant = non_malignant, malignant = malignant), file = "data/processed/seurat_split.rds")
message("Saved split Seurat objects to data/processed/seurat_split.rds")
