#!/usr/bin/env Rscript

# 02-create_seurat_object_and_basic_qc.R
library(Seurat)

inp <- readRDS("data/processed/expr_and_meta.rds")
counts <- inp$counts
meta <- inp$meta

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = "tirosh")

# Add percent.mt (if mt genes exist)
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Basic QC plots and thresholds -- customizable
# Save QC metadata
saveRDS(seurat_obj, file = "data/processed/seurat_before_filtering.rds")
message("Saved Seurat object before filtering: data/processed/seurat_before_filtering.rds")
