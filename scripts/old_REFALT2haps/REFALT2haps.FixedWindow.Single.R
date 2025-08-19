#!/usr/bin/env Rscript

# REFALT2haps Fixed Window - Haplotype Frequency Estimation
# This script estimates haplotype frequencies using a fixed window size
# and outputs both frequency estimates AND distinguishability quality flag
# 
# Usage: Rscript REFALT2haps.FixedWindow.Single.R <chr> <parfile> <mydir> <window_size_kb> <euchrom_start> <euchrom_end>
# Example: Rscript REFALT2haps.FixedWindow.Single.R chr2R helpfiles/JUICE_haplotype_parameters.R process/test 50 5398184 24684540

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  cat("Usage: Rscript REFALT2haps.FixedWindow.Single.R <chr> <parfile> <mydir> <window_size_kb> <euchrom_start> <euchrom_end>\n")
  cat("Example: Rscript REFALT2haps.FixedWindow.Single.R chr2R helpfiles/JUICE_haplotype_parameters.R process/test 50 5398184 24684540\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
window_size_kb <- as.numeric(args[4])
euchrom_start <- as.numeric(args[5])
euchrom_end <- as.numeric(args[6])

cat("=== Fixed Window Haplotype Estimation ===\n")
cat("Chromosome:", mychr, "\n")
cat("Window size:", window_size_kb, "kb\n")
cat("Euchromatin:", euchrom_start, "to", euchrom_end, "bp\n")
cat("Parameter file:", parfile, "\n")
cat("Output directory:", mydir, "\n\n")

# Source unified function
source("scripts/haplotype_estimation_functions.R")

# Load parameters
if (!file.exists(parfile)) {
  cat("Error: Parameter file not found:", parfile, "\n")
  quit(status = 1)
}
source(parfile)

cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), samples (", length(names_in_bam), ")\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")
cat("✓ Samples:", paste(names_in_bam, collapse=", "), "\n\n")

# Load REFALT data
refalt_file <- file.path(mydir, paste0("RefAlt.", mychr, ".txt"))
if (!file.exists(refalt_file)) {
  cat("Error: REFALT file not found:", refalt_file, "\n")
  quit(status = 1)
}

cat("Loading REFALT data from:", refalt_file, "\n")
df <- readr::read_tsv(refalt_file, show_col_types = FALSE)
cat("✓ Raw data:", nrow(df), "rows\n")

# Transform data and calculate frequencies
df2 <- df %>%
  mutate(
    freq = REF / (REF + ALT),
    N = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions and include sample data
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

cat("✓ Quality-filtered positions:", length(quality_filtered_positions), "\n")
cat("✓ Data ready:", nrow(df3), "rows\n")

# Calculate proper grid-based scan positions
# Start: first position >= euchrom_start divisible by step
# End: last position <= euchrom_end divisible by step
scan_start <- ceiling(euchrom_start / step) * step
scan_end <- floor(euchrom_end / step) * step

scan_positions <- seq(scan_start, scan_end, by = step)

cat("✓ Scan grid: ", scan_start, "to", scan_end, "by", step, "\n")
cat("✓ Scan positions:", length(scan_positions), "\n\n")

# Convert window size to bp
window_size_bp <- window_size_kb * 1000

# Professional approach: use expand_grid + purrr::pmap_dfr
cat("Running fixed window haplotype estimation...\n")

results_df <- expand_grid(
  pos = scan_positions,
  sample_name = names_in_bam
) %>%
  purrr::pmap_dfr(~ estimate_haplotypes(
    pos = ..1,
    sample_name = ..2,
    df3 = df3,
    founders = founders,
    h_cutoff = h_cutoff,
    method = "fixed",
    window_size_bp = window_size_bp,
    chr = mychr,
    verbose = 0  # Silent for production
  ))

# Save results
if (nrow(results_df) > 0) {
  # Create output directory if it doesn't exist
  dir.create(mydir, recursive = TRUE, showWarnings = FALSE)
  
  # Save to RDS file
  output_file <- file.path(mydir, paste0("fixed_window_", window_size_kb, "kb_results_", mychr, ".RDS"))
  saveRDS(results_df, output_file)
  
  cat("✓ Results saved:", output_file, "\n")
  cat("✓ Results:", nrow(results_df), "rows\n")
  cat("✓ Positions:", length(unique(results_df$pos)), "unique positions\n")
  cat("✓ Samples:", length(unique(results_df$sample)), "unique samples\n")
  
  # Summary statistics
  success_rate <- mean(results_df$estimate_OK == 1, na.rm = TRUE) * 100
  cat("✓ Success rate:", sprintf("%.1f%%", success_rate), "\n")
  
} else {
  cat("✗ No results generated\n")
  quit(status = 1)
}

cat("\n=== Fixed Window Estimation Complete ===\n")