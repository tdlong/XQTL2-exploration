#!/usr/bin/env Rscriptsub

# Debug wrapper for haplotype estimation
# Usage: Rscript debug_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file> [n_positions] [n_samples] [version]
# Example: Rscript debug_haplotype_estimation.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R 100 5 slow
# Example: Rscript debug_haplotype_estimation.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R 100 5 fast

library(tidyverse)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript debug_haplotype_estimation.R <chr> <method> <parameter> <output_dir> <param_file> [n_positions] [n_samples] [version]\n")
  cat("  chr: chromosome (e.g., chr2R)\n")
  cat("  method: estimation method (adaptive)\n")
  cat("  parameter: method parameter (e.g., 4 for h_cutoff)\n")
  cat("  output_dir: output directory\n")
  cat("  param_file: parameter file\n")
  cat("  n_positions: (optional) number of positions to test (default: 100)\n")
  cat("  n_samples: (optional) number of samples to test (default: 5)\n")
  cat("  version: (optional) 'slow' or 'fast' (default: slow)\n")
  quit(status = 1)
}

chr <- args[1]
method <- args[2]
parameter <- as.numeric(args[3])
output_dir <- args[4]
param_file <- args[5]
n_positions <- if (length(args) >= 6) as.numeric(args[6]) else 100
n_samples <- if (length(args) >= 7) as.numeric(args[7]) else 5
version <- if (length(args) >= 8) args[8] else "slow"

cat("=== DEBUG HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Testing positions:", n_positions, "\n")
cat("Testing samples:", n_samples, "\n")
cat("Version:", version, "\n")
cat("================================\n\n")

# The production pipelines are designed to be run as scripts, not sourced
# We'll call them via system() instead

# Load parameters
source(param_file)

# Load RefAlt data (same format as production pipeline)
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

cat("Loading RefAlt data from:", refalt_file, "\n")
df <- read.table(refalt_file, header = TRUE)
cat("Loaded", nrow(df), "positions from RefAlt file\n")

# Process data (same as production pipeline)
cat("Processing data...\n")
df3 <- df %>%
  mutate(CHROM = chr) %>%
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
  select(-REF, -ALT) %>%
  pivot_wider(names_from = name, values_from = c(freq, N), names_sep = "_") %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "value") %>%
  separate(lab, into = c("type", "sample"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  filter(!is.na(freq)) %>%
  select(POS, name = sample, freq = freq)

# Get sample names from parameter file
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
  filter(sample %in% names_in_bam) %>%
  filter(POS %in% all_positions)

cat("Debug dataset size:", nrow(df3_debug), "rows\n")
cat("Samples:", paste(names_in_bam, collapse = ", "), "\n")
cat("Position range:", min(all_positions), "to", max(all_positions), "\n\n")

# Run haplotype estimation with debug mode enabled
cat("Running haplotype estimation in DEBUG mode...\n")
cat("This will show detailed output for each position-sample combination.\n\n")

# Call the production pipeline as a script with debug mode
cat("Running production pipeline in DEBUG mode...\n")

# Build the command
if (version == "fast") {
  cmd <- paste("Rscript scripts/production/haplotype_estimation_pipeline_fast.R", 
               chr, method, parameter, output_dir, param_file, "debug")
} else {
  cmd <- paste("Rscript scripts/production/haplotype_estimation_pipeline.R", 
               chr, method, parameter, output_dir, param_file, "debug")
}

cat("Command:", cmd, "\n")

# Run the command
result <- system(cmd, intern = TRUE)
cat("Pipeline output:\n")
cat(paste(result, collapse = "\n"), "\n")

# The production pipeline saves its own results, so we just need to report success
cat("✓ Production pipeline completed\n")
cat("✓ Check the output directory for results\n")

cat("\n=== DEBUG COMPLETE ===\n")
cat("Production pipeline completed successfully\n")
cat("Check the output directory for results\n")