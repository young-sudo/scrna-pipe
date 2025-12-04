#!/usr/bin/env Rscript

# 02-create_seurat_object_and_basic_qc.R
library(Seurat)
library(optparse)
library(ggplot2)

inp <- readRDS("results/processed/expr_and_meta.rds")
out_dir <- "results/processed/qc"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=inp, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
inp <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- inp$counts
meta <- inp$meta

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = "tirosh")

# Add percent.mt (if mt genes exist)
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Plot QC figures
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Violin plots for features
vln_plot <- VlnPlot(
  tirosh_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(
  filename = file.path(fig_dir, "qc_violin_plot.png"),
  plot = vln_plot,
  width = 8, height = 4, dpi = 300
)

scatter_plot <- FeatureScatter(
  tirosh_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) + ggtitle(NULL)

ggsave(
  filename = file.path(fig_dir, "count_vs_features.png"),
  plot = scatter_plot,
  width = 5, height = 4, dpi = 300
)


# Basic QC plots and thresholds -- customizable
# Save QC metadata
outfile <- file.path(out_dir, "seurat_qc.rds")
saveRDS(seurat_obj, file = outfile)
message("Saved Seurat object before filtering: ", outfile)
