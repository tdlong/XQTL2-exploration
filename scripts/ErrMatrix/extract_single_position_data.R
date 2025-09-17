#!/usr/bin/env Rscript

# Extract df3 and args for single position - modified from working production code
# This ensures we use the exact same data processing pipeline

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript scripts/ErrMatrix/extract_single_position_data.R <chr> <param_file> <output_dir> <position> <sample_name> <h_cutoff>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
testing_position <- as.numeric(args[4])
sample_name <- args[5]
h_cutoff <- as.numeric(args[6])

cat("=== Extract single position data ===\n")
cat("chr:", chr, "\n")
cat("position:", testing_position, "\n")
cat("sample:", sample_name, "\n")
cat("h_cutoff:", h_cutoff, "\n\n")

# Load parameters
source(param_file)
if (!exists("founders")) stop("Param file must define 'founders'")

# Load the working production functions
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")

# Load RefAlt data (exact same as production)
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

df3 <- process_refalt_data(refalt_file, founders)
cat("✓ df3 loaded:", nrow(df3), "rows x", ncol(df3), "cols\n")

# Pre-subset df3 to df4 for this position (exact same as production)
max_window <- 500000  # Largest window size in est_haps_var
window_start <- max(1, testing_position - max_window/2)
window_end <- testing_position + max_window/2

df4 <- df3 %>%
  dplyr::filter(POS >= window_start & POS <= window_end)

cat("✓ df4 subsetted:", nrow(df4), "SNPs in window", window_start, "-", window_end, "\n")

# Check if sample exists in data
if (!sample_name %in% names(df4)) {
  stop("Sample not found in data: ", sample_name)
}

# Package everything needed for est_haps_var()
payload <- list(
  df3 = df4,  # Use pre-subsetted data like production
  args = list(
    testing_position = testing_position,
    sample_name = sample_name,
    founders = founders,
    h_cutoff = h_cutoff,
    method = "adaptive",
    window_size_bp = NULL,
    chr = chr,
    verbose = 2
  )
)

outfile <- sprintf("single_position_data_%s_%d_%s_h%d.rds", chr, testing_position, sample_name, h_cutoff)
saveRDS(payload, outfile)
cat("Saved:", outfile, "\n")
cat("=== Extraction complete ===\n")
