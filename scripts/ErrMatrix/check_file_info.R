#!/usr/bin/env Rscript

library(tidyverse)

# Check the 4 files
files <- c(
  "process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS",
  "process/ZINC2/haplotype_results_list_format/smooth_h4_results_chr3R.RDS",
  "process/ZINC2/haplotype_results_list_format/smooth_h4/R.haps.chr3R.out.rds",
  "process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds"
)

cat("=== FILE COMPARISON ===\n")

for (file in files) {
  if (file.exists(file)) {
    data <- readRDS(file)
    cat("\n", file, ":\n")
    cat("  Rows:", nrow(data), "\n")
    cat("  First position:", data$pos[1], "\n")
    cat("  Last position:", data$pos[nrow(data)], "\n")
    cat("  Unique positions:", length(unique(data$pos)), "\n")
  } else {
    cat("\n", file, ": FILE NOT FOUND\n")
  }
}