#!/usr/bin/env Rscript

# Test Adaptive Window Distinguishability  
# Tests the core adaptive distinguishability algorithm
# This should match the production script logic exactly

library(tidyverse)

cat("=== ADAPTIVE WINDOW DISTINGUISHABILITY TEST ===\n\n")

# 1. Load parameters (same as production)
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and transform data (same as production)
cat("\n2. Loading and transforming data...\n")
filein <- "process/JUICE/RefAlt.chr2R.txt"
df <- read.table(filein, header = TRUE)

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

# Transform to wide format and apply quality filter ONCE (same as production)
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

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")

# Keep only high-quality positions in the full dataset
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("✓ Data ready:", nrow(df3), "rows\n")

# 3. Test the adaptive window distinguishability
cat("\n3. Testing adaptive window distinguishability...\n")

# Test parameters
test_pos <- 15000000
test_h_cutoff <- 4  # Test specific h_cutoff
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test h_cutoff:", test_h_cutoff, "\n")
cat("✓ Test sample:", test_sample, "(note: sample not used for distinguishability)\n\n")

# Adaptive window expansion: start small, grow until distinguishable or max size
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)  # Progressive sizes
final_window_size <- NA
estimate_OK <- 0
final_n_snps <- 0

cat("Testing window sizes:", paste(window_sizes/1000, "kb", collapse=", "), "\n\n")

for (window_size in window_sizes) {
  cat("--- Testing window size:", window_size/1000, "kb ---\n")
  
  # Define window boundaries
  window_start <- test_pos - window_size/2
  window_end <- test_pos + window_size/2
  
  # Get SNPs in window for founders only (data is already quality-filtered)
  window_snps_long <- df3 %>%
    filter(POS >= window_start & POS <= window_end & name %in% founders)
  
  # Convert to wide format (rows = positions, columns = founders)
  wide_data <- window_snps_long %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  cat("SNPs in window:", nrow(wide_data), "positions\n")
  
  # Check if we have all founder columns and enough data
  if (ncol(wide_data) < length(founders) + 1) {
    cat("✗ Missing founders, trying larger window\n\n")
    next  # Try larger window
  }
  
  if (nrow(wide_data) < 10) {
    cat("✗ Insufficient data in window, trying larger window\n\n")
    next  # Try larger window
  }
  
  # Get founder matrix (no need for additional quality filtering)
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  
  # Convert to matrix for clustering
  founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
  
  if (nrow(founder_matrix_clean) < 10) {
    cat("✗ Insufficient clean data, trying larger window\n\n")
    next  # Try larger window
  }
  
  cat("Founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  
  # Hierarchical clustering to check distinguishability
  tryCatch({
    distances <- dist(t(founder_matrix_clean), method = "euclidean")
    hclust_result <- hclust(distances, method = "ward.D2")
    
    # Cut tree at h_cutoff
    groups <- cutree(hclust_result, h = test_h_cutoff)
    n_groups <- length(unique(groups))
    
    cat("Number of groups at h_cutoff", test_h_cutoff, ":", n_groups, "\n")
    cat("Group assignments:", paste(groups, collapse=", "), "\n")
    
    # Check if all 8 founders can be distinguished
    if (n_groups == 8) {
      cat("✓ SUCCESS: All 8 founders distinguished!\n")
      final_window_size <- window_size
      estimate_OK <- 1
      final_n_snps <- nrow(wide_data)
      break  # Success! Stop expanding
    } else {
      cat("✗ Only", n_groups, "groups, trying larger window\n\n")
    }
    
  }, error = function(e) {
    cat("✗ Clustering failed:", e$message, ", trying larger window\n\n")
    # Clustering failed, try larger window
  })
}

# Record result (either successful at some window size, or failed at all sizes)
if (is.na(final_window_size)) {
  final_window_size <- max(window_sizes)  # Used largest window but failed
  cat("✗ FAILED: Could not distinguish all founders even at", final_window_size/1000, "kb\n")
} else {
  cat("✅ SUCCESS: All founders distinguished at", final_window_size/1000, "kb\n")
}

cat("\n=== FINAL RESULT ===\n")
cat("final_window_size:", final_window_size, "bp (", final_window_size/1000, "kb)\n")
cat("estimate_OK:", estimate_OK, "\n")
cat("n_snps:", final_n_snps, "\n")

cat("\n=== ADAPTIVE WINDOW TEST COMPLETE ===\n")
