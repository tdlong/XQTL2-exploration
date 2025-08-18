#!/usr/bin/env Rscript

# Test Fixed Window Distinguishability
# Tests the core distinguishability algorithm
# This should match the production script logic exactly

library(tidyverse)

cat("=== FIXED WINDOW DISTINGUISHABILITY TEST ===\n\n")

# 1. Load parameters (same as production)
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")
cat("✓ H_cutoff:", h_cutoff, "\n")

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

# 3. Test the fixed window distinguishability
cat("\n3. Testing fixed window distinguishability...\n")

# Test parameters
test_pos <- 15000000
test_window_size <- 50000  # Fixed 50kb window
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test window size:", test_window_size, "bp\n")
cat("✓ Test sample:", test_sample, "(note: sample not used for distinguishability)\n\n")

# Calculate window boundaries
window_start <- test_pos - test_window_size/2
window_end <- test_pos + test_window_size/2
cat("Window range:", window_start, "to", window_end, "bp\n")

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
  cat("✗ Missing founders in window, estimate_OK = 0\n")
  quit()
}

if (nrow(wide_data) < 10) {
  cat("✗ Insufficient data in window, estimate_OK = 0\n")
  quit()
}

# Get founder matrix (no need for additional quality filtering)
founder_matrix <- wide_data %>%
  select(all_of(founders)) %>%
  as.matrix()

# Convert to matrix for clustering
founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]

if (nrow(founder_matrix_clean) < 10) {
  cat("✗ Insufficient clean data for clustering, estimate_OK = 0\n")
  quit()
}

cat("Final founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")

# Hierarchical clustering to check distinguishability
tryCatch({
  distances <- dist(t(founder_matrix_clean), method = "euclidean")
  hclust_result <- hclust(distances, method = "ward.D2")
  
  # Cut tree at h_cutoff
  groups <- cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))
  
  cat("Number of founder groups at h_cutoff", h_cutoff, ":", n_groups, "\n")
  cat("Group assignments:", paste(groups, collapse=", "), "\n")
  
  # Check if all 8 founders can be distinguished
  estimate_OK <- ifelse(n_groups == 8, 1, 0)
  
  cat("\n=== RESULT ===\n")
  cat("estimate_OK:", estimate_OK, "\n")
  
  if (estimate_OK == 1) {
    cat("✓ All 8 founders can be distinguished!\n")
  } else {
    cat("✗ Only", n_groups, "groups - founders cannot be fully distinguished\n")
  }
  
}, error = function(e) {
  cat("✗ Clustering failed:", e$message, "\n")
  cat("estimate_OK: 0\n")
})

cat("\n=== FIXED WINDOW TEST COMPLETE ===\n")
