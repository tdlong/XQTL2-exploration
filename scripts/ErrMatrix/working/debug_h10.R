#!/usr/bin/env Rscript

# Debug h_cutoff=10 to see what causes the catastrophic failure

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# Load the existing hunk data (h_cutoff=4)
hunk_file <- "hunk_data_chr3R_19780000_Rep01_W_F_h4.rds"
if (!file.exists(hunk_file)) {
  stop("Hunk data file not found: ", hunk_file)
}

hunk_data <- readRDS(hunk_file)
cat("Loaded hunk data with", nrow(hunk_data$df3), "rows\n\n")

# Source the est_haps_var function from the debug wrapper
source("scripts/ErrMatrix/working/debug_wrapper.R", local = TRUE)

# Run with h_cutoff=10 and maximum verbosity
cat("=== RUNNING DEBUG WITH H_CUTOFF=10 ===\n")
args_with_h10 <- hunk_data$args
args_with_h10$h_cutoff <- 10
args_with_h10$verbose <- 3  # Maximum verbosity

result <- do.call(est_haps_var, c(list(df3 = hunk_data$df3), args_with_h10))
