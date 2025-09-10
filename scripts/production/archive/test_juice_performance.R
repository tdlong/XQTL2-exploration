#!/usr/bin/env Rscript

# Minimal JUICE performance test
# Usage: Rscript test_juice_performance.R <chr> <n_positions> <n_samples> [version]
# Example: Rscript test_juice_performance.R chr2R 50 3 slow
# Example: Rscript test_juice_performance.R chr2R 50 3 fast

library(tidyverse)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript test_juice_performance.R <chr> <n_positions> <n_samples> [version]\n")
  cat("  chr: chromosome (e.g., chr2R)\n")
  cat("  n_positions: number of positions to test (e.g., 50)\n")
  cat("  n_samples: number of samples to test (e.g., 3)\n")
  cat("  version: 'slow' or 'fast' (default: slow)\n")
  quit(status = 1)
}

chr <- args[1]
n_positions <- as.numeric(args[2])
n_samples <- as.numeric(args[3])
version <- if (length(args) >= 4) args[4] else "slow"

cat("=== JUICE PERFORMANCE TEST ===\n")
cat("Chromosome:", chr, "\n")
cat("Positions:", n_positions, "\n")
cat("Samples:", n_samples, "\n")
cat("Version:", version, "\n")
cat("==============================\n\n")

# Source the appropriate pipeline
if (version == "fast") {
  cat("Loading FAST (vectorized) pipeline...\n")
  source("scripts/production/haplotype_estimation_pipeline_fast.R")
} else {
  cat("Loading SLOW (original) pipeline...\n")
  source("scripts/production/haplotype_estimation_pipeline.R")
}

# Load JUICE parameters
source("helpfiles/JUICE_haplotype_parameters.R")

# Check if JUICE RefAlt data exists
refalt_file <- paste0("process/JUICE/RefAlt.", chr, ".txt")
if (!file.exists(refalt_file)) {
  stop("JUICE RefAlt file not found: ", refalt_file)
}

cat("Found JUICE RefAlt file:", refalt_file, "\n")

# Load a small sample of the data to test
cat("Loading sample data...\n")
df <- read.table(refalt_file, header = TRUE, nrows = n_positions * 10)  # Load more than needed
cat("Loaded", nrow(df), "rows from RefAlt file\n")

# Get first n_positions
positions <- sort(unique(df$POS))[1:n_positions]
cat("Using positions:", min(positions), "to", max(positions), "\n")

# Filter to our test positions
df_test <- df %>%
  filter(POS %in% positions)

# Get first n_samples
samples <- names_in_bam[1:min(n_samples, length(names_in_bam))]
cat("Using samples:", paste(samples, collapse = ", "), "\n")

# Process data (same as production pipeline)
cat("Processing data...\n")
df3 <- df_test %>%
  mutate(CHROM = chr) %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  filter(lab %in% samples) %>%
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
  select(CHROM, POS, sample = sample, freq = freq, N = N)

cat("Processed data:", nrow(df3), "rows\n")

# Quality filter
founder_wide <- df3 %>%
  select(-N) %>%
  pivot_wider(names_from = sample, values_from = freq) %>%
  select(-CHROM, -POS)

# Ensure column order matches founders
founder_wide <- founder_wide[, c("POS", founders)]

# Quality filter: keep only positions where all founders are "fixed"
quality_filter <- apply(founder_wide[, founders], 1, function(row) {
  all(row < 0.03 | row > 0.97)
})

founder_matrix_clean <- founder_wide[quality_filter, ]
cat("Quality-filtered positions:", nrow(founder_matrix_clean), "\n")

# Get positions to scan
scan_positions <- sort(unique(founder_matrix_clean$POS))
cat("Final positions for testing:", length(scan_positions), "\n\n")

# Time the haplotype estimation
cat("Starting haplotype estimation...\n")
start_time <- Sys.time()

if (version == "fast") {
  # Use the vectorized approach
  results <- estimate_haplotypes_vectorized(df3, founders, h_cutoff, debug = FALSE)
} else {
  # Use the original approach with limited data
  results <- expand_grid(
    pos = scan_positions,
    sample_name = samples
  ) %>%
    purrr::pmap_dfr(function(pos, sample_name) {
      result <- estimate_haplotypes_list_format(
        pos = pos,
        sample_name = sample_name,
        df3 = df3,
        founders = founders,
        h_cutoff = h_cutoff,
        window_size_bp = NULL,
        method = "adaptive",
        chr = chr,
        debug = FALSE
      )
      
      return(tibble(
        pos = pos,
        sample = sample_name,
        Groups = list(result$Groups),
        Haps = list(result$Haps),
        Err = list(result$Err),
        Names = list(result$Names)
      ))
    })
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("=== PERFORMANCE RESULTS ===\n")
cat("Version:", version, "\n")
cat("Positions tested:", length(scan_positions), "\n")
cat("Samples tested:", length(samples), "\n")
cat("Total combinations:", length(scan_positions) * length(samples), "\n")
cat("Time elapsed:", round(elapsed, 2), "seconds\n")
cat("Time per combination:", round(elapsed / (length(scan_positions) * length(samples)), 4), "seconds\n")

# Estimate full chromosome time
total_positions <- 2000  # Rough estimate for chr2R
total_samples <- length(names_in_bam)
total_combinations <- total_positions * total_samples
estimated_time <- (elapsed / (length(scan_positions) * length(samples))) * total_combinations

cat("\n=== EXTRAPOLATED ESTIMATES ===\n")
cat("Estimated positions in chr2R:", total_positions, "\n")
cat("Total samples in JUICE:", total_samples, "\n")
cat("Total combinations:", total_combinations, "\n")
cat("Estimated time for full chr2R:", round(estimated_time / 60, 1), "minutes\n")
cat("Estimated time for full chr2R:", round(estimated_time / 3600, 1), "hours\n")

# Save test results
output_file <- paste0("test_juice_", version, "_", chr, "_", n_positions, "pos_", n_samples, "samp.RDS")
saveRDS(results, output_file)
cat("\nTest results saved to:", output_file, "\n")

cat("\n=== TEST COMPLETE ===\n")
