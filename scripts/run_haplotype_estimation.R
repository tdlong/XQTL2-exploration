#!/usr/bin/env Rscript

# Production Haplotype Estimation Wrapper
# Runs chromosome-wide haplotype estimation for a specific method and parameter combination

# Usage: Rscript run_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file>
# Example: Rscript run_haplotype_estimation.R chr2R fixed 50 <output_dir> <param_file>
# Example: Rscript run_haplotype_estimation.R chr2R adaptive 6 <output_dir> <param_file>

library(tidyverse)

# Source unified function
source("scripts/haplotype_estimation_functions.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  cat("Usage: Rscript run_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file>\n")
  cat("Examples:\n")
  cat("  Fixed window:   Rscript run_haplotype_estimation.R chr2R fixed 50 <output_dir> <param_file>\n")
  cat("  Adaptive window: Rscript run_haplotype_estimation.R chr2R adaptive 6 <output_dir> <param_file>\n")
  quit(status = 1)
}

chr <- args[1]
method <- args[2] 
parameter <- as.numeric(args[3])
output_dir <- args[4]
param_file <- args[5]

cat("=== PRODUCTION HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Validate inputs
if (!method %in% c("fixed", "adaptive")) {
  cat("Error: Method must be 'fixed' or 'adaptive'\n")
  quit(status = 1)
}

if (!chr %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")) {
  cat("Error: Chromosome must be one of: chr2L, chr2R, chr3L, chr3R, chrX\n")
  quit(status = 1)
}

# Hardcoded euchromatin boundaries (stolen from Slurm script)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchrom_start <- euchromatin_boundaries[[chr]][1]
euchrom_end <- euchromatin_boundaries[[chr]][2]

cat("Euchromatin boundaries:", euchrom_start, "to", euchrom_end, "bp\n\n")

# 1. Load parameters (same as test wrapper)
cat("1. Loading parameters...\n")

if (!file.exists(param_file)) {
  cat("Error: Parameter file not found:", param_file, "\n")
  quit(status = 1)
}

source(param_file)
cat("✓ Parameter file:", param_file, "\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), samples (", length(names_in_bam), ")\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and transform data (same as test wrapper)
cat("\n2. Loading and transforming data...\n")
filein <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))

if (!file.exists(filein)) {
  cat("Error: REFALT file not found:", filein, "\n")
  quit(status = 1)
}

cat("\nMemory before data load:\n")
print(gc())
df <- read.table(filein, header = TRUE)
cat("\nMemory after data load:\n")
print(gc())

cat("\nMemory before data transform:\n")
print(gc())
df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    N = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Transform to wide format and apply quality filter ONCE (same as test wrapper)
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

# Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions and include sample data
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")
cat("✓ Data ready:", nrow(df3), "rows\n")
cat("\nMemory after data preparation:\n")
print(gc())

# Clear intermediate objects
rm(df2, founder_wide, quality_filtered_positions)
gc()

# 3. Define scan positions (chromosome-wide grid)
cat("\n3. Defining scan positions...\n")

# Calculate proper grid-based scan positions
scan_start <- ceiling(euchrom_start / step) * step
scan_end <- floor(euchrom_end / step) * step
scan_positions <- seq(scan_start, scan_end, by = step)

cat("✓ Scan grid:", scan_start, "to", scan_end, "by", step, "\n")
cat("✓ Scan positions:", length(scan_positions), "\n")

# 4. Run haplotype estimation (same pattern as test wrapper)
cat("\n4. Running haplotype estimation...\n")

# Set up method-specific parameters
if (method == "fixed") {
  window_size_bp <- parameter * 1000  # Convert kb to bp
  h_cutoff_used <- h_cutoff  # Use default from parameter file
  cat("Fixed window size:", parameter, "kb (", window_size_bp, "bp)\n")
} else {
  window_size_bp <- NULL  # Not used for adaptive
  h_cutoff_used <- parameter
  cat("Adaptive h_cutoff:", parameter, "\n")
}

# Use expand_grid + purrr::pmap_dfr (same as test wrapper)
cat("Processing", length(scan_positions), "positions ×", length(names_in_bam), "samples...\n")
cat("\nMemory before processing positions:\n")
print(gc())

results_df <- expand_grid(
  pos = scan_positions,
  sample_name = names_in_bam
) %>%
  purrr::pmap_dfr(~ estimate_haplotypes(
    pos = ..1,
    sample_name = ..2,
    df3 = df3,
    founders = founders,
    h_cutoff = h_cutoff_used,
    method = method,
    window_size_bp = window_size_bp,
    chr = chr,
    verbose = 1  # Show important diagnostics
  ))

# 5. Save results with intelligent filename
cat("\n5. Saving results...\n")

# Generate intelligent output filename
if (method == "fixed") {
  output_filename <- paste0("fixed_window_", parameter, "kb_results_", chr, ".RDS")
} else {
  output_filename <- paste0("adaptive_window_h", parameter, "_results_", chr, ".RDS")
}

# Create results subdirectory
results_dir <- file.path(output_dir, "haplotype_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(results_dir, output_filename)

if (nrow(results_df) > 0) {
  saveRDS(results_df, output_file)
  
  cat("✓ Results saved:", output_file, "\n")
  cat("✓ Results:", nrow(results_df), "rows\n")
  cat("✓ Positions:", length(unique(results_df$pos)), "unique positions\n")
  cat("✓ Samples:", length(unique(results_df$sample)), "unique samples\n")
  
  # Summary statistics
  success_rate <- mean(results_df$estimate_OK == 1, na.rm = TRUE) * 100
  na_rate <- mean(is.na(results_df$estimate_OK)) * 100
  
  cat("✓ Success rate:", sprintf("%.1f%%", success_rate), "\n")
  cat("✓ NA rate:", sprintf("%.1f%%", na_rate), "\n")
  
} else {
  cat("✗ No results generated\n")
  quit(status = 1)
}

cat("\n=== PRODUCTION HAPLOTYPE ESTIMATION COMPLETE ===\n")
