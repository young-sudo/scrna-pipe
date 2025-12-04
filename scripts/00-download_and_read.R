#!/usr/bin/env Rscript

# 00-download_and_read.R
# Download (optional) and read the raw matrix file into R
# Configure paths

library(optparse)

# raw_dir <- "data/raw"
out_dir <- "results/processed/read"
file_in <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=file_in, help="Input file"),
  make_option(c("-o", "--outdir"), type="character", default=out_dir, help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

file_in <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file_in)) {
  message("File not found: ", file_in)
  message("Please download from GEO:")
  message("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz")
  stop("Data file missing")
}

# Read table (fast) -- stringsAsFactors = FALSE
message("Reading data...")
raw <- data.table::fread(file_in, header = TRUE, data.table = FALSE)

file_out <- file.path(out_dir, "raw_table.rds")
# Save a copy in RDS for faster loading later
saveRDS(raw, file = file_out)
message("Saved raw RDS to ", file_out)
