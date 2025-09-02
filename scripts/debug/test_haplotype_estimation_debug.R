#!/usr/bin/env Rscript

# Test Haplotype Estimation Debug
# This script tests the haplotype estimation process for a specific position
# without modifying production code
#
# Usage: Rscript test_haplotype_estimation_debug.R <chr> <param_file> <output_dir> <estimator> <position>

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript test_haplotype_estimation_debug.R <chr> <param_file> <output_dir> <estimator> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]
target_position <- as.numeric(args[5])

cat("=== Test Haplotype Estimation Debug ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Estimator:", estimator, "\n")
cat("Target position:", target_position, "\n\n")

# Load parameter file
cat("Loading parameter file...\n")
source(param_file)
cat("✓ Parameter file loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Samples to process:", paste(names_in_bam, collapse = ", "), "\n\n")

# Define euchromatin boundaries
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]

# Load REFALT data
cat("Loading REFALT data...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT file not found:", refalt_file)
}

# Load and process REFALT data
df <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies
cat("Converting counts to frequencies...\n")
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

cat("✓ Processed REFALT data:", nrow(df2), "rows\n")

# Filter for high-quality SNPs
cat("Filtering for high-quality SNPs...\n")
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for evaluation (euchromatin only)
valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n\n")

# Now test the haplotype estimation process for the target position
cat("=== Testing Haplotype Estimation at Position", target_position, "===\n")

# Determine window size for fixed method
if (grepl("^fixed_", estimator)) {
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  window_size_bp <- window_size_kb * 1000
  cat("Fixed window size:", window_size_kb, "kb (", window_size_bp, "bp)\n")
} else {
  cat("Adaptive method - using 20kb for testing\n")
  window_size_bp <- 20000
}

# Calculate window boundaries
window_start <- target_position - window_size_bp/2
window_end <- target_position + window_size_bp/2

cat("Window boundaries:", window_start, "-", window_end, "bp\n")
cat("Window center:", target_position, "bp\n\n")

# Get data in the window
cat("Getting data in window...\n")
window_data <- valid_snps %>%
  filter(POS >= window_start & POS <= window_end & name %in% c(founders, names_in_bam[1]))

cat("Total SNPs in window:", nrow(window_data %>% distinct(CHROM, POS)), "\n")

if (nrow(window_data) == 0) {
  cat("❌ No data in window - cannot proceed\n")
  quit(status = 1)
}

# Create wide format data
wide_data <- window_data %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("Wide data dimensions:", nrow(wide_data), "x", ncol(wide_data), "\n")
cat("Columns:", paste(names(wide_data), collapse = ", "), "\n\n")

# Check if we have all required data
if (!all(c(founders, names_in_bam[1]) %in% names(wide_data))) {
  missing_cols <- setdiff(c(founders, names_in_bam[1]), names(wide_data))
  cat("❌ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  quit(status = 1)
}

if (nrow(wide_data) < 10) {
  cat("❌ Insufficient SNPs in window:", nrow(wide_data), "< 10\n")
  quit(status = 1)
}

# Get founder matrix and sample frequencies
founder_matrix <- wide_data %>%
  select(all_of(founders)) %>%
  as.matrix()

sample_freqs <- wide_data %>%
  pull(!!names_in_bam[1])

cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
cat("Sample frequencies length:", length(sample_freqs), "\n\n")

# Remove incomplete rows
complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
sample_freqs_clean <- sample_freqs[complete_rows]

cat("Complete SNPs for analysis:", nrow(founder_matrix_clean), "\n")

if (nrow(founder_matrix_clean) < 10) {
  cat("❌ Too few complete SNPs:", nrow(founder_matrix_clean), "< 10\n")
  quit(status = 1)
}

# Now show detailed debugging information
cat("\n=== DETAILED DEBUGGING BEFORE LSEI ===\n")
cat("Founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
cat("Sample frequencies length:", length(sample_freqs_clean), "\n")
cat("h_cutoff: 2.5 (default)\n")

# Show founder matrix summary
cat("\nFounder matrix summary:\n")
for (i in seq_along(founders)) {
  founder_data <- founder_matrix_clean[, i]
  cat(sprintf("  %s: %d SNPs, range %.3f-%.3f, mean %.3f\n", 
             founders[i], length(founder_data), 
             min(founder_data, na.rm=TRUE), max(founder_data, na.rm=TRUE),
             mean(founder_data, na.rm=TRUE)))
}

# Show sample frequency summary
cat("\nSample frequency summary:\n")
cat("  Range:", min(sample_freqs_clean, na.rm=TRUE), "-", max(sample_freqs_clean, na.rm=TRUE), "\n")
cat("  Mean:", mean(sample_freqs_clean, na.rm=TRUE), "\n")
cat("  NAs:", sum(is.na(sample_freqs_clean)), "\n")

# Calculate and show founder distances
cat("\nFounder distance matrix:\n")
founder_dist <- dist(t(founder_matrix_clean))
founder_dist_matrix <- as.matrix(founder_dist)
print(round(founder_dist_matrix, 4))

# Show clustering results
cat("\nHierarchical clustering results:\n")
h_cutoff <- 2.5  # Default value
hclust_result <- hclust(founder_dist, method = "complete")
groups <- cutree(hclust_result, h = h_cutoff)
n_groups <- length(unique(groups))

cat("  h_cutoff:", h_cutoff, "\n")
cat("  Number of groups:", n_groups, "\n")
cat("  Expected groups:", length(founders), "\n")
cat("  Groups sufficient:", n_groups == length(founders), "\n")

# Show group assignments
unique_clusters <- unique(groups)
for (cluster_id in unique_clusters) {
  cluster_founders <- founders[groups == cluster_id]
  if (length(cluster_founders) == 1) {
    cat(sprintf("  Group %d: %s (individual)\n", cluster_id, cluster_founders))
  } else {
    cat(sprintf("  Group %d: %s\n", cluster_id, paste(cluster_founders, collapse=", ")))
  }
}

# Show minimum distance between any two founders
min_dist <- min(founder_dist)
cat("\nDistance analysis:\n")
cat("  Minimum distance between any two founders:", round(min_dist, 4), "\n")
cat("  h_cutoff satisfied:", min_dist >= h_cutoff, "\n")

# Show which founders are closest
min_dist_idx <- which(founder_dist_matrix == min_dist, arr.ind = TRUE)
if (nrow(min_dist_idx) > 0) {
  closest_pair <- min_dist_idx[1, ]
  cat("  Closest founders:", founders[closest_pair[1]], "and", founders[closest_pair[2]], "\n")
}

cat("=== END DEBUGGING ===\n\n")

# Now try to run LSEI
cat("=== Attempting LSEI ===\n")
tryCatch({
  E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
  F <- 1.0
  G <- diag(length(founders))  # Non-negativity constraints
  H <- matrix(rep(0.0003, length(founders)))  # Lower bound
  
  cat("LSEI inputs:\n")
  cat("Matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  cat("Number of constraints:", nrow(E), "equality +", nrow(G), "inequality\n")
  
  # Try LSEI
  lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                               E = E, F = F, G = G, H = H)
  
  if (lsei_result$IsError == 0) {
    cat("✓ LSEI successful!\n")
    cat("Haplotype frequency estimates:\n")
    haplotype_freqs <- lsei_result$X
    names(haplotype_freqs) <- founders
    for (i in seq_along(founders)) {
      cat(sprintf("  %s: %.4f\n", founders[i], haplotype_freqs[i]))
    }
    cat("Sum of frequencies:", round(sum(haplotype_freqs), 4), "\n")
  } else {
    cat("✗ LSEI failed with error code:", lsei_result$IsError, "\n")
  }
  
}, error = function(e) {
  cat("✗ LSEI error:", e$message, "\n")
})

cat("\n=== Test Complete ===\n")
cat("This shows exactly what data went into the haplotype estimation.\n")
cat("Check if the window boundaries and clustering criteria are correct.\n")
