#!/usr/bin/env Rscript

# Extract specific sample data for direct comparison with debug result
# Target: Rep01_W_F at position 19780000

suppressPackageStartupMessages({
  library(tidyverse)
})

data_file <- "testing_positions_comparison.rds"

if (!file.exists(data_file)) {
  stop("Data file not found: ", normalizePath(data_file))
}

cat("Loading data from:", data_file, "\n")
xx <- readRDS(data_file)

# Unnest the parallel list-columns to per-sample rows
per_sample <- xx %>%
  dplyr::select(CHROM, pos, method, sample, Haps, Err) %>%
  tidyr::unnest(c(sample, Haps, Err)) %>%
  dplyr::rename(sample_id = sample, hap_vec = Haps, err_mat = Err)

# Filter for our target position and sample
target_data <- per_sample %>%
  dplyr::filter(pos == 19780000 & sample_id == "Rep01_W_F")

if (nrow(target_data) == 0) {
  stop("No data found for Rep01_W_F at position 19780000")
}

cat("Found", nrow(target_data), "rows for Rep01_W_F at position 19780000\n\n")

# Extract and display the data
for (i in 1:nrow(target_data)) {
  method <- target_data$method[i]
  haps <- as.numeric(target_data$hap_vec[[i]])
  err_mat <- as.matrix(target_data$err_mat[[i]])
  err_diag_sum <- sum(diag(err_mat))
  
  cat("=== METHOD:", method, "===\n")
  cat("Haplotype frequencies:\n")
  cat(paste(sprintf("%.6f", haps), collapse = ", "), "\n")
  cat("Error matrix diagonal sum:", err_diag_sum, "\n")
  cat("Error matrix diagonal values:\n")
  cat(paste(sprintf("%.2e", diag(err_mat)), collapse = ", "), "\n\n")
}

# Also save to file for easy comparison
output_data <- target_data %>%
  dplyr::mutate(
    hap_frequencies = purrr::map_chr(hap_vec, ~paste(sprintf("%.6f", as.numeric(.x)), collapse = ", ")),
    err_diag_sum = purrr::map_dbl(err_mat, ~sum(diag(as.matrix(.x)))),
    err_diag_values = purrr::map_chr(err_mat, ~paste(sprintf("%.2e", diag(as.matrix(.x))), collapse = ", "))
  ) %>%
  dplyr::select(method, hap_frequencies, err_diag_sum, err_diag_values)

readr::write_csv(output_data, "Rep01_W_F_19780000_comparison.csv")
cat("Saved detailed comparison to: Rep01_W_F_19780000_comparison.csv\n")
