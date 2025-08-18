#!/usr/bin/env Rscript

# Test Fixed Window Distinguishability
# Tests the core distinguishability algorithm
# This should match the production script logic exactly

library(tidyverse)
library(limSolve)

cat("=== FIXED WINDOW DISTINGUISHABILITY TEST ===\n\n")

# 1. Load parameters (same as production)
cat("1. Loading parameters...\n")
source("helpfiles/JUICE_haplotype_parameters.R")
cat("✓ Parameter file: helpfiles/JUICE_haplotype_parameters.R\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), h_cutoff (", h_cutoff, "), samples (", length(names_in_bam), ")\n")
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
# Use samples defined in parameter file
non_founder_samples <- names_in_bam

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

# Get SNPs in window for both founders AND sample (data is already quality-filtered)
window_snps_long <- df3 %>%
  filter(POS >= window_start & POS <= window_end & name %in% c(founders, test_sample))

# Convert to wide format (rows = positions, columns = founders + sample)
wide_data <- window_snps_long %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("SNPs in window:", nrow(wide_data), "positions\n")
cat("Columns available:", paste(names(wide_data), collapse = ", "), "\n")

# Check if we have all founder columns + sample and enough data
if (!all(c(founders, test_sample) %in% names(wide_data))) {
  missing_cols <- c(founders, test_sample)[!c(founders, test_sample) %in% names(wide_data)]
  cat("✗ Missing columns in window:", paste(missing_cols, collapse = ", "), "\n")
  quit()
}

if (nrow(wide_data) < 10) {
  cat("✗ Insufficient data in window, estimate_OK = 0\n")
  quit()
}

# Get founder matrix and sample frequencies for LSEI
founder_matrix <- wide_data %>%
  select(all_of(founders)) %>%
  as.matrix()

sample_freqs <- wide_data %>%
  pull(!!test_sample)

# Remove rows with any NAs (both founder and sample data must be complete)
complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
sample_freqs_clean <- sample_freqs[complete_rows]

if (nrow(founder_matrix_clean) < 10) {
  cat("✗ Insufficient clean data for clustering, estimate_OK = 0\n")
  quit()
}

cat("Final founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")

# Show some raw data for diagnostics
cat("\n=== FOUNDER FREQUENCY RANGES ===\n")
for (i in 1:ncol(founder_matrix_clean)) {
  founder_name <- founders[i]
  freq_range <- range(founder_matrix_clean[, i], na.rm = TRUE)
  cat(founder_name, ": ", sprintf("%.3f - %.3f", freq_range[1], freq_range[2]), "\n")
}

# Show sample of raw SNP data
cat("\n=== SAMPLE SNP DATA (first 10 positions, percentages) ===\n")
sample_data <- wide_data[1:min(10, nrow(wide_data)), ]

# Create formatted table
cat(sprintf("%-10s", "POS"), paste(sprintf("%-3s", founders), collapse=" "), "\n")
cat(paste(rep("-", 10 + length(founders) * 4), collapse=""), "\n")

for (i in 1:nrow(sample_data)) {
  pos <- sample_data$POS[i]
  freqs <- sprintf("%2.0f", as.numeric(sample_data[i, founders]) * 100)
  cat(sprintf("%-10s", pos), paste(sprintf("%3s", freqs), collapse=" "), "\n")
}

# Initialize estimate_OK variable
estimate_OK <- NA

# Hierarchical clustering to check distinguishability
cat("\n=== CLUSTERING ANALYSIS ===\n")
tryCatch({
  distances <- dist(t(founder_matrix_clean), method = "euclidean")
  hclust_result <- hclust(distances, method = "ward.D2")
  
  # Show clustering distances as a matrix
  cat("Clustering distance matrix:\n")
  dist_matrix <- as.matrix(distances)
  rownames(dist_matrix) <- founders
  colnames(dist_matrix) <- founders
  
  # Print matrix header
  cat(sprintf("%4s", ""), paste(sprintf("%6s", founders), collapse=""), "\n")
  
  # Print matrix rows
  for (i in 1:length(founders)) {
    row_values <- sprintf("%6.1f", dist_matrix[i, ])
    cat(sprintf("%4s", founders[i]), paste(row_values, collapse=""), "\n")
  }
  
  # Cut tree at h_cutoff
  groups <- cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))
  
  cat("\nClustering results at h_cutoff", h_cutoff, ":\n")
  cat("Number of groups:", n_groups, "\n")
  cat("Group assignments:", paste(groups, collapse=", "), "\n")
  
  # Show which founders are in which groups
  for (group_id in unique(groups)) {
    group_founders <- founders[groups == group_id]
    cat("Group", group_id, ":", paste(group_founders, collapse=", "), "\n")
  }
  
  # Run LSEI first to get actual haplotype frequency estimates
  cat("\n=== LSEI HAPLOTYPE ESTIMATION ===\n")
  tryCatch({
    # LSEI constraints: sum to 1, non-negative
    E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    G <- diag(length(founders))  # Non-negativity constraints
    H <- matrix(rep(0.0003, length(founders)))  # Lower bound
    
    lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                 E = E, F = F, G = G, H = H)
    
    if (lsei_result$IsError == 0) {
      cat("✓ LSEI successful!\n")
      haplotype_freqs <- lsei_result$X
      names(haplotype_freqs) <- founders
      
      cat("Haplotype frequency estimates:\n")
      for (i in seq_along(founders)) {
        cat(sprintf("  %s: %.4f\n", founders[i], haplotype_freqs[i]))
      }
      cat("Sum of frequencies:", round(sum(haplotype_freqs), 4), "\n")
      
      # Now check distinguishability
      estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
      
      # Show which founders are in which groups with their frequencies
      cat("\nGroup-wise frequencies:\n")
      for (group_id in unique(groups)) {
        group_founders <- founders[groups == group_id]
        group_freq <- sum(haplotype_freqs[groups == group_id])
        cat("Group", group_id, ":", paste(group_founders, collapse="+"), "=", round(group_freq, 4), "\n")
      }
    } else {
      cat("✗ LSEI failed with error code:", lsei_result$IsError, "\n")
      estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
    }
  }, error = function(e) {
    cat("✗ LSEI error:", e$message, "\n")
    estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
  })
  
  cat("\n=== FINAL RESULT ===\n")
  cat("estimate_OK:", estimate_OK, "\n")
  
  if (estimate_OK == 1) {
    cat("✓ All", length(founders), "founders can be distinguished individually!\n")
    cat("✓ Perfect distinguishability at", test_window_size/1000, "kb window size\n")
    cat("✓ LSEI estimation successful - haplotype frequencies available\n")
  } else if (estimate_OK == 0) {
    cat("✗ Only", n_groups, "groups - some founders cannot be distinguished\n")
    cat("✗ Need larger window or different h_cutoff for full distinguishability\n")
    cat("! LSEI may still provide estimates, but they may be unreliable\n")
  } else {
    cat("✗ LSEI estimation failed - no haplotype frequencies available\n")
  }
  
}, error = function(e) {
  cat("✗ Clustering failed:", e$message, "\n")
  estimate_OK <- 0
  cat("estimate_OK:", estimate_OK, "\n")
})

cat("\n=== FIXED WINDOW TEST COMPLETE ===\n")
