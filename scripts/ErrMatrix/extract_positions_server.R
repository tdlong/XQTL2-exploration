#!/usr/bin/env Rscript

# Simple extraction script for server execution
# Run this on the server to extract the testing positions

library(tidyverse)

# Define the testing positions
testing_positions <- c(19780000, 19790000, 19800000, 19810000, 19820000, 19830000, 19840000)

# File paths
adapt_file <- "/dfs7/adl/tdlong/exploration/XQTL2-exploration/process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds"
fixed_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/R.haps.chr3R.out.rds"

# Extract adaptive data
cat("Loading adaptive data...\n")
adapt_data <- readRDS(adapt_file)
if (is.list(adapt_data) && !is.data.frame(adapt_data)) {
  # Find the data frame with 'pos' column
  for (i in seq_along(adapt_data)) {
    if (is.data.frame(adapt_data[[i]]) && "pos" %in% names(adapt_data[[i]])) {
      adapt_data <- adapt_data[[i]]
      break
    }
  }
}
adapt_extracted <- adapt_data %>%
  dplyr::filter(pos %in% testing_positions) %>%
  dplyr::mutate(method = "adapt")

# Extract fixed data
cat("Loading fixed data...\n")
fixed_data <- readRDS(fixed_file)
if (is.list(fixed_data) && !is.data.frame(fixed_data)) {
  # Find the data frame with 'pos' column
  for (i in seq_along(fixed_data)) {
    if (is.data.frame(fixed_data[[i]]) && "pos" %in% names(fixed_data[[i]])) {
      fixed_data <- fixed_data[[i]]
      break
    }
  }
}
fixed_extracted <- fixed_data %>%
  dplyr::filter(pos %in% testing_positions) %>%
  dplyr::mutate(method = "fixed")

# Combine and save
combined <- dplyr::bind_rows(adapt_extracted, fixed_extracted) %>%
  dplyr::arrange(pos, method)

saveRDS(combined, "testing_positions_comparison.rds")

cat("Extraction complete!\n")
cat("Found", nrow(adapt_extracted), "adaptive positions\n")
cat("Found", nrow(fixed_extracted), "fixed positions\n")
cat("Saved to: testing_positions_comparison.rds\n")
