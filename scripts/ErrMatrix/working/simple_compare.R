#!/usr/bin/env Rscript

# Simple comparison of hunk result vs production data

suppressPackageStartupMessages({
  library(tidyverse)
})

# Load our debug result
debug_result <- readRDS("debug_result_chr3R_19780000_Rep01_W_F.rds")
cat("=== DEBUG RESULT ===\n")
cat("Error matrix diagonal sum:", sum(diag(debug_result$Err)), "\n")
cat("Haplotype estimates:", paste(sprintf("%.6f", debug_result$Haps), collapse = ", "), "\n\n")

# Load production data
prod_data <- readRDS("testing_positions_comparison.rds")
cat("=== PRODUCTION DATA ===\n")
cat("Total rows:", nrow(prod_data), "\n")

# Find Rep01_W_F at position 19780000
pos_data <- prod_data %>%
  dplyr::filter(pos == 19780000 & sample == "Rep01_W_F")

if (nrow(pos_data) == 0) {
  cat("No data found for Rep01_W_F at position 19780000\n")
  quit()
}

cat("Found", nrow(pos_data), "rows for Rep01_W_F at 19780000\n")

for (i in 1:nrow(pos_data)) {
  method <- pos_data$method[i]
  err_sum <- sum(diag(pos_data$Err[[i]]))
  haps <- pos_data$Haps[[i]]
  
  cat("\nMethod:", method, "\n")
  cat("Error diagonal sum:", err_sum, "\n")
  cat("Haplotype estimates:", paste(sprintf("%.6f", haps), collapse = ", "), "\n")
  
  if (method == "adapt") {
    cat("\n=== COMPARISON: DEBUG vs PRODUCTION ADAPT ===\n")
    cat("Error sum difference (debug - prod):", sum(diag(debug_result$Err)) - err_sum, "\n")
    cat("Error sum ratio (debug / prod):", sum(diag(debug_result$Err)) / err_sum, "\n")
    
    hap_diff <- sum(abs(debug_result$Haps - haps))
    cat("Haplotype difference (sum of abs diffs):", hap_diff, "\n")
  }
}
