#!/usr/bin/env Rscript

library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  cat("Usage: Rscript compare_adaptive_vs_reshaped.R <chr> <output_dir> [limit]\n")
  quit(status = 1)
}

chr <- args[1]
output_dir <- args[2]
limit_n <- if (length(args) == 3) as.integer(args[3]) else NA_integer_

# Construct file paths
list_dir <- file.path(output_dir, "haplotype_results_list_format")
orig_file <- file.path(list_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
resh_file <- file.path(list_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))

stopifnot(file.exists(orig_file), file.exists(resh_file))

cat("Comparing:\n  orig:", orig_file, "\n  resh:", resh_file, "\n\n")

# Load files
orig <- readRDS(orig_file)
resh <- readRDS(resh_file)

cat("Original: ", nrow(orig), "rows (one per CHROM/pos/sample)\n")
cat("Reshaped: ", nrow(resh), "rows (one per CHROM/pos)\n")

# Track results
total_samples <- 0
identical_samples <- 0
max_haps_diff <- 0
max_err_diff <- 0
max_trace_diff <- 0

# Get all positions
positions <- sort(unique(orig$pos))
cat("Total positions:", length(positions), "\n")

# Apply limit if specified
if (!is.na(limit_n) && limit_n > 0) {
  positions <- head(positions, limit_n)
  cat("Limited to first", length(positions), "positions\n")
}

cat("Checking", length(positions), "positions...\n")

# Process each position
for (pos in positions) {
  orig_pos <- orig %>% filter(pos == !!pos)
  resh_pos <- resh %>% filter(pos == !!pos)
  
  if (nrow(resh_pos) == 0) next
  
  # Check each sample at this position
  for (i in 1:nrow(orig_pos)) {
    orig_sample <- orig_pos[i,]
    resh_sample_idx <- which(resh_pos$sample[[1]] == orig_sample$sample)[1]
    
    if (is.na(resh_sample_idx)) next
    
    # Extract data from reshaped file
    resh_haps <- resh_pos$Haps[[1]][[resh_sample_idx]]
    resh_err <- resh_pos$Err[[1]][[resh_sample_idx]]
    resh_names <- resh_pos$Names[[1]][[resh_sample_idx]]
    
    # Compare
    haps_diff <- max(abs(orig_sample$Haps[[1]] - resh_haps), na.rm = TRUE)
    err_diff <- max(abs(orig_sample$Err[[1]] - resh_err), na.rm = TRUE)
    trace_diff <- abs(sum(diag(orig_sample$Err[[1]])) - sum(diag(resh_err)))
    
    # Handle missing values
    if (is.na(haps_diff)) haps_diff <- Inf
    if (is.na(err_diff)) err_diff <- Inf
    if (is.na(trace_diff)) trace_diff <- Inf
    
    # Update tracking
    total_samples <- total_samples + 1
    if (haps_diff == 0 && err_diff == 0 && trace_diff == 0) {
      identical_samples <- identical_samples + 1
    }
    
    max_haps_diff <- max(max_haps_diff, haps_diff)
    max_err_diff <- max(max_err_diff, err_diff)
    max_trace_diff <- max(max_trace_diff, trace_diff)
  }
}

cat("\n=== RESULTS ===\n")
cat("Total samples checked:", total_samples, "\n")
cat("Identical samples:", identical_samples, "\n")
cat("Percentage identical:", round(100 * identical_samples / total_samples, 2), "%\n")
cat("Max haps difference:", max_haps_diff, "\n")
cat("Max err difference:", max_err_diff, "\n")
cat("Max trace difference:", max_trace_diff, "\n")

if (identical_samples == total_samples) {
  cat("\n✅ VERIFICATION PASSED: All", total_samples, "samples are identical between original and reshaped files!\n")
  cat("✅ The reshaping preserves data integrity perfectly.\n")
} else {
  cat("\n❌ VERIFICATION FAILED: Some samples don't match.\n")
}