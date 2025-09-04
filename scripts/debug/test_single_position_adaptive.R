#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# =============================================================================
# Test Adaptive Algorithm on Single Position
# =============================================================================
# 
# This script tests the adaptive algorithm on a single position to diagnose
# why it might be failing to produce valid results.
#
# USAGE:
# Rscript scripts/debug/test_single_position_adaptive.R <chr> <param_file> <output_dir> <position>
#
# EXAMPLE:
# Rscript scripts/debug/test_single_position_adaptive.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 5400000
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/debug/test_single_position_adaptive.R <chr> <param_file> <output_dir> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])

cat("=== TESTING ADAPTIVE ALGORITHM ON SINGLE POSITION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Load observed data
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found:", refalt_file)
}

cat("Loading REFALT data from:", refalt_file, "\n")
refalt_data <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                       col_types = "cdcdd", show_col_types = FALSE)

cat("✓ REFALT data loaded:", nrow(refalt_data), "rows\n")

# Convert counts to frequencies
observed_data <- refalt_data %>%
  mutate(
    total_count = ref_count + alt_count,
    freq = alt_count / total_count
  ) %>%
  filter(total_count > 0) %>%
  select(chr, pos, name, freq)

cat("✓ Converted to frequencies:", nrow(observed_data), "rows\n")

# Filter for euchromatin
euchromatin_boundaries <- list(
  chr2L = c(1, 23011544),
  chr2R = c(5398184, 24684540),
  chr3L = c(1, 24543557),
  chr3R = c(1, 27905053),
  chr4 = c(1, 1351857),
  chrX = c(1, 22422827)
)

if (chr %in% names(euchromatin_boundaries)) {
  boundaries <- euchromatin_boundaries[[chr]]
  observed_data <- observed_data %>%
    filter(pos >= boundaries[1] & pos <= boundaries[2])
  cat("✓ Filtered for euchromatin:", nrow(observed_data), "rows\n")
  cat("Euchromatin region:", boundaries[1], "-", boundaries[2], "bp\n")
}

# Get samples
samples <- unique(observed_data$name)
samples <- samples[!samples %in% founders]
cat("Samples found:", paste(samples, collapse = ", "), "\n\n")

# Test on first sample
test_sample <- samples[1]
cat("Testing on sample:", test_sample, "\n\n")

# Test the adaptive algorithm with different window sizes
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)

for (window_size in window_sizes) {
  cat("=== TESTING WINDOW SIZE:", window_size, "bp ===\n")
  
  window_start <- test_position - window_size/2
  window_end <- test_position + window_size/2
  
  window_data <- observed_data %>%
    filter(pos >= window_start & pos <= window_end & name %in% c(founders, test_sample))
  
  cat("SNPs in window:", nrow(window_data), "\n")
  
  if (nrow(window_data) == 0) {
    cat("❌ No SNPs found in window\n\n")
    next
  }
  
  # Convert to wide format
  wide_data <- window_data %>%
    select(pos, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  if (!all(c(founders, test_sample) %in% names(wide_data))) {
    cat("❌ Missing founders or sample in wide data\n\n")
    next
  }
  
  # Get founder matrix and sample frequencies
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  sample_freqs <- wide_data %>%
    pull(!!test_sample)
  
  # Remove rows with missing data
  complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  cat("Clean data dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  
  if (nrow(founder_matrix_clean) < 10) {
    cat("❌ Insufficient SNPs after cleaning:", nrow(founder_matrix_clean), "\n\n")
    next
  }
  
  # Test hierarchical clustering
  cat("Testing hierarchical clustering...\n")
  founder_dist <- dist(t(founder_matrix_clean))
  hclust_result <- hclust(founder_dist, method = "complete")
  groups <- cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))
  
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  cat("Number of groups:", n_groups, "\n")
  
  # Test LSEI
  cat("Testing LSEI...\n")
  n_founders <- ncol(founder_matrix_clean)
  E <- matrix(rep(1, n_founders), nrow = 1)
  F <- 1.0
  
  result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                          E = E, F = F, 
                          G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                          fulloutput = TRUE)
  
  cat("LSEI IsError:", result$IsError, "\n")
  
  if (result$IsError == 0) {
    cat("✓ LSEI successful\n")
    hap_freqs <- result$X
    names(hap_freqs) <- founders
    cat("Haplotype frequencies sum:", round(sum(hap_freqs), 4), "\n")
    cat("Error matrix available:", !is.null(result$cov), "\n")
    
    if (n_groups == 8) {
      cat("✓ All 8 founders distinguishable!\n")
    } else {
      cat("⚠️  Only", n_groups, "groups found\n")
    }
  } else {
    cat("❌ LSEI failed with error code:", result$IsError, "\n")
  }
  
  cat("\n")
}

cat("=== TEST COMPLETE ===\n")
