#!/usr/bin/env Rscript

# Debug wrapper for haplotype estimation
# Usage: Rscript debug_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file> [n_positions] [n_samples]
# Example: Rscript debug_haplotype_estimation.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R 100 5

library(tidyverse)

# Source the main pipeline
source("scripts/production/haplotype_estimation_pipeline.R")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript debug_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file> [n_positions] [n_samples]\n")
  cat("  chr: chromosome (e.g., chr2R)\n")
  cat("  method: estimation method (adaptive)\n")
  cat("  parameter: method parameter (e.g., 4 for h_cutoff)\n")
  cat("  output_dir: output directory\n")
  cat("  param_file: parameter file\n")
  cat("  n_positions: (optional) number of positions to test (default: 100)\n")
  cat("  n_samples: (optional) number of samples to test (default: 5)\n")
  quit(status = 1)
}

chr <- args[1]
method <- args[2]
parameter <- as.numeric(args[3])
output_dir <- args[4]
param_file <- args[5]
n_positions <- if (length(args) >= 6) as.numeric(args[6]) else 100
n_samples <- if (length(args) >= 7) as.numeric(args[7]) else 5

cat("=== DEBUG HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Testing positions:", n_positions, "\n")
cat("Testing samples:", n_samples, "\n")
cat("================================\n\n")

# Load parameters
source(param_file)

# Load RefAlt data
refalt_file <- file.path(output_dir, paste0("RefAlt_", chr, ".RDS"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

cat("Loading RefAlt data from:", refalt_file, "\n")
df3 <- readRDS(refalt_file)

# Get sample names
names_in_bam <- unique(df3$lab)
cat("Total samples available:", length(names_in_bam), "\n")

# Limit to first n_samples for debugging
if (length(names_in_bam) > n_samples) {
  names_in_bam <- names_in_bam[1:n_samples]
  cat("Limited to first", n_samples, "samples for debugging\n")
}

# Get positions
all_positions <- unique(df3$POS)
cat("Total positions available:", length(all_positions), "\n")

# Limit to first n_positions for debugging
if (length(all_positions) > n_positions) {
  all_positions <- all_positions[1:n_positions]
  cat("Limited to first", n_positions, "positions for debugging\n")
}

# Filter data for debugging
df3_debug <- df3 %>%
  filter(lab %in% names_in_bam) %>%
  filter(POS %in% all_positions)

cat("Debug dataset size:", nrow(df3_debug), "rows\n")
cat("Samples:", paste(names_in_bam, collapse = ", "), "\n")
cat("Position range:", min(all_positions), "to", max(all_positions), "\n\n")

# Run haplotype estimation with debug mode enabled
cat("Running haplotype estimation in DEBUG mode...\n")
cat("This will show detailed output for each position-sample combination.\n\n")

# Create a modified version of the function that uses our debug dataset
run_haplotype_estimation_debug <- function(chr, method, parameter, output_dir, param_file) {
  # Load parameters
  source(param_file)
  
  # Use our debug dataset instead of loading from file
  df3 <- df3_debug
  
  # Get founders from the data
  founders <- unique(df3$founder)
  founders <- founders[!is.na(founders)]
  
  cat("Founders:", paste(founders, collapse = ", "), "\n")
  cat("Number of founders:", length(founders), "\n\n")
  
  # Quality filter
  founder_wide <- df3 %>%
    select(POS, founder, lab, freq) %>%
    pivot_wider(names_from = founder, values_from = freq, values_fill = 0) %>%
    filter(lab %in% names_in_bam)
  
  quality_filter <- apply(founder_wide[, founders], 1, function(row) {
    all(row < 0.03 | row > 0.97)
  })
  
  founder_matrix_clean <- founder_wide[quality_filter, ]
  cat("Quality-filtered positions:", nrow(founder_matrix_clean), "\n")
  
  # Get positions to scan
  scan_positions <- sort(unique(founder_matrix_clean$POS))
  cat("Using", length(scan_positions), "positions for haplotype estimation\n\n")
  
  # Run haplotype estimation for each position and sample
  cat("Running haplotype estimation...\n")
  total_combinations <- length(scan_positions) * length(names_in_bam)
  cat("Processing", total_combinations, "position-sample combinations\n\n")
  
  results <- expand_grid(
    pos = scan_positions,
    sample_name = names_in_bam
  ) %>%
    mutate(combo_id = row_number()) %>%
    purrr::pmap_dfr(function(pos, sample_name, combo_id) {
      cat("=== Position", pos, "Sample", sample_name, "===\n")
      
      result <- estimate_haplotypes_list_format(
        pos = pos,
        sample_name = sample_name,
        df3 = df3,
        founders = founders,
        h_cutoff = parameter,
        window_size_bp = NULL,
        method = method,
        chr = chr,
        debug = TRUE  # Enable debug mode
      )
      
      cat("Result: Groups =", paste(result$Groups, collapse = " "), "\n")
      cat("Result: Haps =", paste(round(result$Haps, 3), collapse = " "), "\n")
      cat("\n")
      
      return(tibble(
        pos = pos,
        sample = sample_name,
        Groups = list(result$Groups),
        Haps = list(result$Haps),
        Err = list(result$Err),
        Names = list(result$Names)
      ))
    })
  
  # Save results
  output_file <- file.path(output_dir, paste0("debug_", chr, "_", method, "_", parameter, ".RDS"))
  cat("Saving debug results to:", output_file, "\n")
  saveRDS(results, output_file)
  
  cat("✓ Debug results saved\n")
  return(results)
}

# Run the debug version
results <- run_haplotype_estimation_debug(chr, method, parameter, output_dir, param_file)

cat("\n=== DEBUG COMPLETE ===\n")
cat("Results saved to:", file.path(output_dir, paste0("debug_", chr, "_", method, "_", parameter, ".RDS")), "\n")
cat("Total combinations processed:", nrow(results), "\n")
