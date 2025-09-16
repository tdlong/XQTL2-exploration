#!/usr/bin/env Rscript

# Debug script to check matrix names in reshaped data
# Usage: Rscript scripts/ErrMatrix/debug_matrix_names.R <chr> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript scripts/ErrMatrix/debug_matrix_names.R <chr> <output_dir>\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]

list_dir <- file.path(output_dir, "haplotype_results_list_format")
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

cat("Loading reshaped file...\n")
resh <- readRDS(resh_file)

# Check first row before unnesting
first_row <- resh[1,]
cat("Before unnesting - first row:\n")
cat("  Err length:", length(first_row$Err[[1]]), "\n")
if (length(first_row$Err[[1]]) > 0) {
  first_err <- first_row$Err[[1]][[1]]
  cat("  First Err is.matrix:", is.matrix(first_err), "\n")
  if (is.matrix(first_err)) {
    cat("  First Err dim:", nrow(first_err), "x", ncol(first_err), "\n")
    cat("  First Err rownames:", paste(rownames(first_err), collapse=","), "\n")
    cat("  First Err colnames:", paste(colnames(first_err), collapse=","), "\n")
  }
}

# Unnest and check
resh_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

cat("\nAfter unnesting - first row:\n")
first_row_u <- resh_u[1,]
err_matrix <- first_row_u$Err[[1]]
cat("  Err is.matrix:", is.matrix(err_matrix), "\n")
if (is.matrix(err_matrix)) {
  cat("  Err dim:", nrow(err_matrix), "x", ncol(err_matrix), "\n")
  cat("  Err rownames:", paste(rownames(err_matrix), collapse=","), "\n")
  cat("  Err colnames:", paste(colnames(err_matrix), collapse=","), "\n")
  cat("  Err rownames is.null:", is.null(rownames(err_matrix)), "\n")
  cat("  Err colnames is.null:", is.null(colnames(err_matrix)), "\n")
}

# Check if the matrix has any attributes
cat("\nMatrix attributes:\n")
print(attributes(err_matrix))

cat("\nDone.\n")
