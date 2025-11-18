#!/usr/bin/env Rscript

# 06-malignant_preprocess_and_dimreduce.R
library(Seurat)

objs <- readRDS("data/processed/seurat_split.rds")
malignant <- objs$malignant

malignant <- NormalizeData(malignant)
all.genes <- rownames(malignant)
malignant <- ScaleData(malignant, features = all.genes)
malignant <- RunPCA(malignant, features = all.genes)

# tSNE/UMAP
set.seed(40)
malignant <- RunTSNE(malignant, dims = 1:15, seed.use = 40)
set.seed(1)
malignant <- RunUMAP(malignant, dims = 1:15, seed.use = 1)

saveRDS(malignant, file = "data/processed/malignant_processed.rds")
message("Saved malignant processed object to data/processed/malignant_processed.rds")
