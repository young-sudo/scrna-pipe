#!/usr/bin/env Rscript

# 03-filter_and_split_malignant.R
library(Seurat)
library(optparse)

seurat_obj <- readRDS("results/processed/seurat_qc.rds")
out_dir <- "results/"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=seurat_obj, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
seurat_obj <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# apply threshold-based filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 20)

# The dataset contains a metadata column:
# malignant.1.no.2.yes.0.unresolved.
# Unfortunately, this is not standard
if (!"malignant.1.no.2.yes.0.unresolved." %in% colnames(seurat_obj@meta.data)) {
  stop("Expected metadata column malignant.1.no.2.yes.0.unresolved. not found")
}

non_malignant <- subset(seurat_obj, subset = malignant.1.no.2.yes.0.unresolved. == 1)
malignant <- subset(seurat_obj, subset = malignant.1.no.2.yes.0.unresolved. == 2)

# saveRDS(list(non_malignant = non_malignant, malignant = malignant), file = outfile)
# Non-malignant
nonmal_dir <- file.path(out_dir, "non_malignant", "seurat")
if (!dir.exists(nonmal_dir)) dir.create(nonmal_dir, recursive = TRUE, showWarnings = TRUE)
saveRDS(non_malignant, file = file.path(nonmal_dir, "non_malignant_seurat.rds"))

# Malignant
mal_dir <- file.path(out_dir, "malignant", "seurat")
if (!dir.exists(mal_dir)) dir.create(mal_dir, recursive = TRUE, showWarnings = TRUE)
saveRDS(malignant, file = file.path(mal_dir, "malignant_seurat.rds"))

message("Saved split Seurat objects to ", out_dir)
