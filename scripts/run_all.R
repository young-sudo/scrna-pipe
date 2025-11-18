#!/usr/bin/env Rscript

# run_all.R

# Master script to run the entire Tirosh melanoma single-cell reanalysis pipeline

# Executed in order from data download to trajectory analysis

# Make sure you are in the repo root

scripts <- list(
"scripts/00-download_and_read.R",
"scripts/01-prepare_expression_and_metadata.R",
"scripts/02-create_seurat_object_and_basic_qc.R",
"scripts/03-filter_and_split_malignant.R",
"scripts/04-nonmalignant_preprocess_and_dimreduce.R",
"scripts/05-nonmalignant_clustering_and_qc_plots.R",
"scripts/06-malignant_preprocess_and_dimreduce.R",
"scripts/07-trajectory_monocle3.R"
)

for (script in scripts) {
cat(sprintf("\n---\nRunning: %s\n", script))
source(script)
cat(sprintf("Finished: %s\n---\n", script))
}

cat("All scripts executed successfully. Check 'data/processed/' and 'figures/' for outputs.\n")
