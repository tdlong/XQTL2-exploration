#!/usr/bin/env Rscript

# Debug script to check unnesting structure
# Usage: Rscript scripts/ErrMatrix/debug_unnesting.R <chr> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript scripts/ErrMatrix/debug_unnesting.R <chr> <output_dir>\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]

list_dir <- file.path(output_dir, "haplotype_results_list_format")
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

cat("Loading reshaped file...\n")
resh <- readRDS(resh_file)

cat("Before unnesting:\n")
cat("  nrow:", nrow(resh), "\n")
cat("  columns:", paste(names(resh), collapse=", "), "\n")

# Check first row structure
first_row <- resh[1,]
cat("\nFirst row before unnesting:\n")
cat("  sample length:", length(first_row$sample[[1]]), "\n")
cat("  Names length:", length(first_row$Names[[1]]), "\n")
cat("  Names content:", paste(first_row$Names[[1]], collapse=","), "\n")
cat("  Err length:", length(first_row$Err[[1]]), "\n")
if (length(first_row$Err[[1]]) > 0) {
  cat("  First Err is.matrix:", is.matrix(first_row$Err[[1]][[1]]), "\n")
  if (is.matrix(first_row$Err[[1]][[1]])) {
    cat("  First Err dim:", nrow(first_row$Err[[1]][[1]]), "x", ncol(first_row$Err[[1]][[1]]), "\n")
    cat("  First Err rownames:", paste(rownames(first_row$Err[[1]][[1]]), collapse=","), "\n")
  }
}

# Unnest
cat("\nUnnesting...\n")
resh_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

cat("After unnesting:\n")
cat("  nrow:", nrow(resh_u), "\n")
cat("  columns:", paste(names(resh_u), collapse=", "), "\n")

# Check first row after unnesting
first_row_u <- resh_u[1,]
cat("\nFirst row after unnesting:\n")
cat("  sample:", first_row_u$sample, "\n")
cat("  Names length:", length(first_row_u$Names[[1]]), "\n")
cat("  Names content:", paste(first_row_u$Names[[1]], collapse=","), "\n")
cat("  Err is.matrix:", is.matrix(first_row_u$Err[[1]]), "\n")
if (is.matrix(first_row_u$Err[[1]])) {
  cat("  Err dim:", nrow(first_row_u$Err[[1]]), "x", ncol(first_row_u$Err[[1]]), "\n")
  cat("  Err rownames:", paste(rownames(first_row_u$Err[[1]]), collapse=","), "\n")
}

cat("\nDone.\n")
