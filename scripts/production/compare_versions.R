#!/usr/bin/env Rscript

# Compare results between slow and fast versions
# Usage: Rscript compare_versions.R <chr> <method> <parameter> <output_dir>
# Example: Rscript compare_versions.R chr2R adaptive 4 process/ZINC2

library(tidyverse)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("Usage: Rscript compare_versions.R <chr> <method> <parameter> <output_dir>\n")
  cat("  chr: chromosome (e.g., chr2R)\n")
  cat("  method: estimation method (adaptive)\n")
  cat("  parameter: method parameter (e.g., 4 for h_cutoff)\n")
  cat("  output_dir: output directory\n")
  quit(status = 1)
}

chr <- args[1]
method <- args[2]
parameter <- as.numeric(args[3])
output_dir <- args[4]

cat("=== COMPARING SLOW vs FAST VERSIONS ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load results from both versions
slow_file <- file.path(output_dir, paste0("debug_slow_", chr, "_", method, "_", parameter, ".RDS"))
fast_file <- file.path(output_dir, paste0("debug_fast_", chr, "_", method, "_", parameter, ".RDS"))

if (!file.exists(slow_file)) {
  stop("Slow version results not found: ", slow_file)
}

if (!file.exists(fast_file)) {
  stop("Fast version results not found: ", fast_file)
}

cat("Loading results...\n")
slow_results <- readRDS(slow_file)
fast_results <- readRDS(fast_file)

cat("Slow version:", nrow(slow_results), "combinations\n")
cat("Fast version:", nrow(fast_results), "combinations\n\n")

# Compare basic structure
cat("=== STRUCTURE COMPARISON ===\n")
cat("Slow version columns:", paste(names(slow_results), collapse = ", "), "\n")
cat("Fast version columns:", paste(names(fast_results), collapse = ", "), "\n")

# Check if they have the same combinations
slow_combos <- slow_results %>%
  select(pos, sample) %>%
  arrange(pos, sample)

fast_combos <- fast_results %>%
  select(pos, sample) %>%
  arrange(pos, sample)

if (identical(slow_combos, fast_combos)) {
  cat("✓ Both versions have identical position-sample combinations\n")
} else {
  cat("✗ Different position-sample combinations!\n")
  cat("Slow combinations:", nrow(slow_combos), "\n")
  cat("Fast combinations:", nrow(fast_combos), "\n")
}

# Compare haplotype frequencies
cat("\n=== HAPLOTYPE FREQUENCY COMPARISON ===\n")

# Extract haplotype frequencies for comparison
slow_haps <- slow_results %>%
  select(pos, sample, Haps) %>%
  mutate(version = "slow")

fast_haps <- fast_results %>%
  select(pos, sample, Haps) %>%
  mutate(version = "fast")

# Combine and compare
combined_haps <- bind_rows(slow_haps, fast_haps) %>%
  group_by(pos, sample) %>%
  summarise(
    slow_haps = list(Haps[version == "slow"][[1]]),
    fast_haps = list(Haps[version == "fast"][[1]]),
    .groups = "drop"
  ) %>%
  mutate(
    hap_diff = map2_dbl(slow_haps, fast_haps, ~ {
      if (length(.x) != length(.y)) return(NA)
      max(abs(.x - .y), na.rm = TRUE)
    })
  )

cat("Maximum haplotype frequency difference:", max(combined_haps$hap_diff, na.rm = TRUE), "\n")
cat("Mean haplotype frequency difference:", mean(combined_haps$hap_diff, na.rm = TRUE), "\n")

# Show some examples
cat("\nFirst 5 position-sample combinations:\n")
print(combined_haps %>%
  select(pos, sample, hap_diff) %>%
  head(5))

# Check for exact matches
exact_matches <- sum(combined_haps$hap_diff < 1e-10, na.rm = TRUE)
total_combinations <- nrow(combined_haps)
cat("\nExact matches (diff < 1e-10):", exact_matches, "/", total_combinations, 
    "(", round(exact_matches/total_combinations*100, 1), "%)\n")

# Check for close matches (diff < 0.01)
close_matches <- sum(combined_haps$hap_diff < 0.01, na.rm = TRUE)
cat("Close matches (diff < 0.01):", close_matches, "/", total_combinations, 
    "(", round(close_matches/total_combinations*100, 1), "%)\n")

# Check for large differences
large_diffs <- sum(combined_haps$hap_diff > 0.1, na.rm = TRUE)
cat("Large differences (diff > 0.1):", large_diffs, "/", total_combinations, 
    "(", round(large_diffs/total_combinations*100, 1), "%)\n")

if (large_diffs > 0) {
  cat("\nPosition-sample combinations with large differences:\n")
  print(combined_haps %>%
    filter(hap_diff > 0.1) %>%
    select(pos, sample, hap_diff) %>%
    arrange(desc(hap_diff)))
}

cat("\n=== COMPARISON COMPLETE ===\n")
