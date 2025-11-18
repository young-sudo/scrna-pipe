#!/usr/bin/env Rscript

# 04-nonmalignant_preprocess_and_dimreduce.R
library(Seurat)

objs <- readRDS("data/processed/seurat_split.rds")
non_malignant <- objs$non_malignant

# Annotate cell types from the metadata mapping used in the Rmd
cluster_names <- c("unresolved", "T-cells", "B-cells", "Macrophages", "Endothelial", "CAF", "NK")
if ("non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK" %in% colnames(non_malignant@meta.data)) {
  non_malignant@meta.data$cell_type <- cluster_names[non_malignant@meta.data$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK + 1]
}

# Normalization / scaling / PCA
non_malignant <- NormalizeData(non_malignant)
all.genes <- rownames(non_malignant)
non_malignant <- ScaleData(non_malignant, features = all.genes)
non_malignant <- RunPCA(non_malignant, features = all.genes)

# Variable features + PCA using variable features
non_malignant <- FindVariableFeatures(non_malignant, selection.method = "vst", nfeatures = 1000)
non_malignant <- ScaleData(non_malignant)
non_malignant <- RunPCA(non_malignant, features = VariableFeatures(object = non_malignant))

# Run UMAP/tSNE using first 15 PCs (set seeds for reproducibility)
set.seed(3)
non_malignant <- RunTSNE(non_malignant, dims = 1:15, seed.use = 3)
set.seed(1)
non_malignant <- RunUMAP(non_malignant, dims = 1:15, seed.use = 1, n.neighbors = 10)

saveRDS(non_malignant, file = "data/processed/non_malignant_processed.rds")
message("Saved non-malignant processed object to data/processed/non_malignant_processed.rds")
