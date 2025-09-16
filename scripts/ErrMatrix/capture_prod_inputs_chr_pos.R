#!/usr/bin/env Rscript

# Capture the actual df4 window and simple summaries that would be passed to
# estimate_haplotypes_list_format inside BASE_VAR_WIDE.R for a single
# chromosome/position/sample, WITHOUT running the full workflow and WITHOUT
# writing any production outputs.
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

# Source production code (functions only; no writes are triggered here)
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Load parameters
source(param_file, local = TRUE)

# Paths to write capture (local, non-production)
capture_rds <- paste0("prod_capture_", chr, "_", target_pos, "_", target_sample, ".RDS")
capture_txt <- paste0("prod_capture_", chr, "_", target_pos, "_", target_sample, ".txt")

# DO NOT run the full workflow. Instead, emulate the exact df4 pre-subset
# used by run_adaptive_estimation and capture that input, safely and quickly.

# Load RefAlt and process -> df3 (wide), identical to BASE_VAR_WIDE.R
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
stopifnot(file.exists(refalt_file))
df3 <- process_refalt_data(refalt_file, founders)

# Pre-subset df3 to df4 (largest 500kb window around target_pos), identical logic
max_window <- 500000
window_start <- max(1, target_pos - max_window/2)
window_end <- target_pos + max_window/2
df4 <- df3 %>% dplyr::filter(POS >= window_start & POS <= window_end)

# Summaries and capture
df4_cols <- names(df4)
df4_n <- nrow(df4)
pos_min <- min(df4$POS, na.rm = TRUE)
pos_max <- max(df4$POS, na.rm = TRUE)
pos_mean <- mean(df4$POS, na.rm = TRUE)
df4_sha1 <- digest(df4, algo = "sha1")

sink(capture_txt)
cat("df4_rows:", df4_n, "\n")
cat("df4_cols:", paste(df4_cols, collapse=","), "\n")
cat("pos_min:", pos_min, " pos_max:", pos_max, " pos_mean:", pos_mean, "\n")
cat("df4_sha1:", df4_sha1, "\n")
cat("founders:", paste(founders, collapse=","), "\n")
cat("h_cutoff:", get("h_cutoff"), " method:", "adaptive", " chr:", chr, "\n")
sink()

saveRDS(list(
  df4 = df4,
  founders = founders,
  h_cutoff = get("h_cutoff"),
  method = "adaptive",
  chr = chr,
  testing_position = target_pos,
  sample_name = target_sample,
  df4_sha1 = df4_sha1
), capture_rds)

cat("Done. Wrote:", capture_txt, "and", capture_rds, "\n")


