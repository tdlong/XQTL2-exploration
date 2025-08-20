#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/find_difficult_regions.R <chr> <param_file> <output_dir>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]

# Load parameters
source(param_file)
results_dir <- file.path(output_dir, "haplotype_results")

cat("=== FINDING DIFFICULT REGIONS ===\n")
cat("Chromosome:", chr, "\n")
cat("Looking for regions where small fixed windows struggle...\n\n")

# Load haplotype results for fixed_20kb (smallest window)
fixed_20kb_file <- file.path(results_dir, "fixed_window_20kb_results_chr2R.RDS")

if (!file.exists(fixed_20kb_file)) {
  stop("Fixed 20kb results not found: ", fixed_20kb_file)
}

results <- readRDS(fixed_20kb_file)

# Calculate estimate_OK rate by position
position_stats <- results %>%
  group_by(pos) %>%
  summarise(
    total_estimates = n(),
    ok_estimates = sum(estimate_OK == 1, na.rm = TRUE),
    failed_estimates = sum(estimate_OK == 0, na.rm = TRUE),
    na_estimates = sum(is.na(estimate_OK)),
    ok_rate = ok_estimates / total_estimates,
    .groups = "drop"
  ) %>%
  arrange(pos)

cat("Overall statistics:\n")
cat("Total positions:", nrow(position_stats), "\n")
cat("Average OK rate:", round(mean(position_stats$ok_rate), 3), "\n")
cat("Positions with OK rate < 0.5:", sum(position_stats$ok_rate < 0.5), "\n")
cat("Positions with OK rate < 0.1:", sum(position_stats$ok_rate < 0.1), "\n\n")

# Find regions with low OK rates (difficult regions)
difficult_regions <- position_stats %>%
  filter(ok_rate < 0.5) %>%
  mutate(
    region_start = floor(pos / 100000) * 100000,  # 100kb regions
    region_end = region_start + 100000
  ) %>%
  group_by(region_start, region_end) %>%
  summarise(
    n_positions = n(),
    avg_ok_rate = mean(ok_rate),
    min_ok_rate = min(ok_rate),
    max_ok_rate = max(ok_rate),
    .groups = "drop"
  ) %>%
  arrange(avg_ok_rate) %>%
  head(10)  # Top 10 most difficult regions

cat("Top 10 most difficult regions (100kb windows):\n")
print(difficult_regions)

# Find specific positions with very low OK rates
worst_positions <- position_stats %>%
  filter(ok_rate < 0.3) %>%
  arrange(ok_rate) %>%
  head(20)

cat("\nTop 20 worst individual positions:\n")
print(worst_positions)

# Save results
output_file <- file.path(results_dir, paste0("difficult_regions_", chr, ".RDS"))
saveRDS(list(
  position_stats = position_stats,
  difficult_regions = difficult_regions,
  worst_positions = worst_positions
), output_file)

cat("\n✓ Results saved to:", output_file, "\n")

# Suggest a region to plot
if (nrow(difficult_regions) > 0) {
  best_region <- difficult_regions[1, ]
  region_center <- (best_region$region_start + best_region$region_end) / 2
  region_center_10kb <- region_center / 10000
  
  cat("\n=== SUGGESTED PLOTTING REGION ===\n")
  cat("Region:", format(best_region$region_start, big.mark=","), "-", format(best_region$region_end, big.mark=","), "\n")
  cat("Center position (10kb units):", region_center_10kb, "\n")
  cat("Average OK rate:", round(best_region$avg_ok_rate, 3), "\n")
  cat("Number of positions:", best_region$n_positions, "\n")
  
  cat("\nCommand to plot this region:\n")
  cat("Rscript scripts/plot_B1_frequencies.R chr2R helpfiles/JUICE_haplotype_parameters.R process/JUICE AJ_1_1\n")
  cat("(Then manually adjust the zoom region in the script to center on", region_center_10kb, ")\n")
}
