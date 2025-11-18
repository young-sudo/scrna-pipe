#!/usr/bin/env Rscript

# 08-utils.R
# Put reusable helper functions here (plot saving, gene checks, etc.)

save_png <- function(plot_expr, filename, width = 800, height = 800) {
  png(filename, width = width, height = height)
  print(plot_expr)
  dev.off()
}

# small helper to check column presence
check_cols <- function(obj, cols) {
  missing <- setdiff(cols, colnames(obj@meta.data))
  if (length(missing) > 0) stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}
