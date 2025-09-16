#!/usr/bin/env Rscript

# Capture the actual df3 (df4 window) and simple summaries passed to
# estimate_haplotypes_list_format inside BASE_VAR_WIDE.R for a single
# chromosome/position/sample, without modifying production code.
#
# Usage:
# Rscript scripts/ErrMatrix/capture_prod_inputs_chr_pos.R <chr> <position> <sample_name> <output_dir> <param_file>

suppressPackageStartupMessages({
  library(tidyverse)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  cat("Usage: Rscript scripts/ErrMatrix/capture_prod_inputs_chr_pos.R <chr> <position> <sample_name> <output_dir> <param_file>\n")
  quit(status = 1)
}

chr <- args[1]
target_pos <- as.numeric(args[2])
target_sample <- args[3]
output_dir <- args[4]
param_file <- args[5]

cat("=== CAPTURE PRODUCTION INPUTS ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", target_pos, "\n")
cat("Sample:", target_sample, "\n")
cat("Output dir:", output_dir, "\n")
cat("Param file:", param_file, "\n\n")

# Source production code
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Load parameters
source(param_file, local = TRUE)

# Keep a reference to original function
original_estimator <- estimate_haplotypes_list_format

# Path to write capture
capture_rds <- paste0("prod_capture_", chr, "_", target_pos, "_", target_sample, ".RDS")
capture_txt <- paste0("prod_capture_", chr, "_", target_pos, "_", target_sample, ".txt")

# Wrap estimator to capture its inputs when pos/sample match
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                                            method = "adaptive", window_size_bp = NULL,
                                            chr = "chr2R", verbose = 0) {
  if (!exists(".captured_flag", inherits = FALSE)) .captured_flag <<- FALSE
  if (!.captured_flag && pos == target_pos && sample_name == target_sample) {
    # Summaries for df3 (this is the pre-subsetted window from run_adaptive_estimation)
    df3_cols <- names(df3)
    df3_n <- nrow(df3)
    pos_min <- min(df3$POS, na.rm = TRUE)
    pos_max <- max(df3$POS, na.rm = TRUE)
    pos_mean <- mean(df3$POS, na.rm = TRUE)
    df3_sha1 <- digest(df3, algo = "sha1")

    cat("[CAPTURE] Matched target; writing capture files\n")
    # Save a compact text summary
    sink(capture_txt)
    cat("df3_rows:", df3_n, "\n")
    cat("df3_cols:", paste(df3_cols, collapse=","), "\n")
    cat("pos_min:", pos_min, " pos_max:", pos_max, " pos_mean:", pos_mean, "\n")
    cat("df3_sha1:", df3_sha1, "\n")
    cat("founders:", paste(founders, collapse=","), "\n")
    cat("h_cutoff:", h_cutoff, " method:", method, " chr:", chr, "\n")
    sink()

    # Save the full object to RDS (df3 + params)
    saveRDS(list(
      df3 = df3,
      founders = founders,
      h_cutoff = h_cutoff,
      method = method,
      chr = chr,
      position = pos,
      sample_name = sample_name,
      df3_sha1 = df3_sha1
    ), capture_rds)

    .captured_flag <<- TRUE
    cat("[CAPTURE] Wrote:", capture_txt, "and", capture_rds, "\n")
  }
  # Delegate to original
  original_estimator(pos, sample_name, df3, founders, h_cutoff, method, window_size_bp, chr, verbose)
}

# Run only the target chromosome in production mode (non-verbose)
cat("Running run_adaptive_estimation for target chromosome...\n")
run_adaptive_estimation(chr = chr, method = "adaptive", parameter = get("h_cutoff"),
                        output_dir = output_dir, param_file = param_file,
                        debug = FALSE, verbose = FALSE, debug_level = 0)

cat("\nDone. If capture matched, see:", capture_txt, "and", capture_rds, "\n")


