#!/usr/bin/env Rscript

# Test Adaptive Window Distinguishability  
# Tests the core adaptive distinguishability algorithm
# This should match the production script logic exactly

library(tidyverse)
library(limSolve)

cat("=== ADAPTIVE WINDOW DISTINGUISHABILITY TEST ===\n\n")

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

# 3. Test the adaptive window distinguishability
cat("\n3. Testing adaptive window distinguishability...\n")

# Test parameters
test_pos <- 15000000
test_h_cutoff <- 10  # Test higher h_cutoff to force window expansion
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test h_cutoff:", test_h_cutoff, "\n")
cat("✓ Test sample:", test_sample, "(note: sample not used for distinguishability)\n\n")

# Adaptive window expansion: start small, grow until distinguishable or max size
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)  # Progressive sizes
final_window_size <- NA
estimate_OK <- 0
final_n_snps <- 0

cat("=== STARTING ADAPTIVE WINDOW ALGORITHM ===\n")
cat("Window sizes to test:", paste(window_sizes/1000, "kb", collapse=", "), "\n\n")

# Track SNP counts across all window sizes
snps_tracking <- data.frame(
  window_size = window_sizes,
  total_snps = 0,
  founder_complete = 0,
  clustering_valid = 0,
  n_groups = 0,
  estimate_OK = 0
)

for (window_idx in seq_along(window_sizes)) {
  window_size <- window_sizes[window_idx]
  cat("--- Window", window_idx, ":", window_size/1000, "kb ---\n")
  
  # Define window boundaries
  window_start <- test_pos - window_size/2
  window_end <- test_pos + window_size/2
  cat("Window range:", window_start, "to", window_end, "bp\n")
  
  # Get SNPs in window for founders only (data is already quality-filtered)
  window_snps_long <- df3 %>%
    filter(POS >= window_start & POS <= window_end & name %in% founders)
  
  # Convert to wide format (rows = positions, columns = founders)
  wide_data <- window_snps_long %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  cat("SNPs in window:", nrow(wide_data), "positions\n")
  
  # Track total SNPs
  snps_tracking$total_snps[window_idx] <- nrow(wide_data)
  
  # Show first 10 SNPs as a table for the first window
  if (window_idx == 1) {
    cat("\n=== FIRST 10 SNPS IN WINDOW (RAW DATA, percentages) ===\n")
    display_data <- wide_data %>%
      arrange(POS) %>%
      head(10)
    
    # Create formatted table
    cat(sprintf("%-10s", "POS"), paste(sprintf("%-3s", founders), collapse=" "), "\n")
    cat(paste(rep("-", 10 + length(founders) * 4), collapse=""), "\n")
    
    for (i in 1:nrow(display_data)) {
      pos <- display_data$POS[i]
      freqs <- sprintf("%2.0f", as.numeric(display_data[i, founders]) * 100)
      cat(sprintf("%-10s", pos), paste(sprintf("%3s", freqs), collapse=" "), "\n")
    }
    cat("\n")
  }
  
  # Check if we have all founder columns and enough data
  if (ncol(wide_data) < length(founders) + 1) {
    cat("✗ Missing founders, trying larger window\n\n")
    snps_tracking$founder_complete[window_idx] <- 0
    next  # Try larger window
  }
  
  if (nrow(wide_data) < 10) {
    cat("✗ Insufficient data in window, trying larger window\n\n")
    snps_tracking$founder_complete[window_idx] <- 0
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
    snps_tracking$founder_complete[window_idx] <- nrow(founder_matrix_clean)
    snps_tracking$clustering_valid[window_idx] <- 0
    next  # Try larger window
  }
  
  snps_tracking$founder_complete[window_idx] <- nrow(founder_matrix)
  snps_tracking$clustering_valid[window_idx] <- nrow(founder_matrix_clean)
  
  cat("Founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  
  # Show founder frequency ranges for this window
  cat("Founder frequency ranges:\n")
  for (i in 1:ncol(founder_matrix_clean)) {
    founder_name <- founders[i]
    freq_range <- range(founder_matrix_clean[, i], na.rm = TRUE)
    cat("  ", founder_name, ": ", sprintf("%.3f - %.3f", freq_range[1], freq_range[2]), "\n")
  }
  
  # Hierarchical clustering to check distinguishability
  tryCatch({
    distances <- dist(t(founder_matrix_clean), method = "euclidean")
    hclust_result <- hclust(distances, method = "ward.D2")
    
    # Show clustering distance matrix
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
    
    min_dist <- min(dist_matrix[dist_matrix > 0])
    max_dist <- max(dist_matrix)
    cat("Min:", sprintf("%.1f", min_dist), "| Max:", sprintf("%.1f", max_dist), "| H_cutoff:", test_h_cutoff, "\n")
    
    # Cut tree at h_cutoff
    groups <- cutree(hclust_result, h = test_h_cutoff)
    n_groups <- length(unique(groups))
    
    snps_tracking$n_groups[window_idx] <- n_groups
    
    cat("Clustering results:\n")
    cat("  Number of groups:", n_groups, "\n")
    cat("  Group assignments:", paste(groups, collapse=", "), "\n")
    
    # Show which founders are in which groups
    for (group_id in unique(groups)) {
      group_founders <- founders[groups == group_id]
      cat("  Group", group_id, ":", paste(group_founders, collapse=", "), "\n")
    }
    
    # Check if all founders can be distinguished
    if (n_groups == length(founders)) {
      cat("✓ SUCCESS: All", length(founders), "founders distinguished individually!\n")
      final_window_size <- window_size
      estimate_OK <- 1
      final_n_snps <- nrow(wide_data)
      snps_tracking$estimate_OK[window_idx] <- 1
      break  # Success! Stop expanding
    } else {
      cat("✗ Only", n_groups, "groups - need larger window\n")
      cat("✗ Some founders still clustered together\n")
      snps_tracking$estimate_OK[window_idx] <- 0
    }
    
  }, error = function(e) {
    cat("✗ Clustering failed:", e$message, ", trying larger window\n")
    snps_tracking$n_groups[window_idx] <- 0
    snps_tracking$estimate_OK[window_idx] <- 0
  })
  
  cat("\n")
}

# Record result (either successful at some window size, or failed at all sizes)
if (is.na(final_window_size)) {
  final_window_size <- max(window_sizes)  # Used largest window but failed
  cat("✗ FAILED: Could not distinguish all founders even at", final_window_size/1000, "kb\n")
} else {
  cat("✅ SUCCESS: All founders distinguished at", final_window_size/1000, "kb\n")
}

# Show comprehensive SNP tracking summary
cat("\n=== SNP TRACKING SUMMARY ===\n")
print(snps_tracking)

# Run LSEI estimation using the optimal window size
cat("\n=== LSEI HAPLOTYPE ESTIMATION ===\n")
if (estimate_OK == 1) {
  # Use the final window size to get data for LSEI
  window_start <- test_pos - final_window_size/2
  window_end <- test_pos + final_window_size/2
  
  window_data <- df3 %>%
    filter(POS >= window_start & POS <= window_end & name %in% c(founders, test_sample))
  
  wide_data <- window_data %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  # Extract founder and sample data
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  sample_freqs <- wide_data %>%
    pull(!!test_sample)
  
  # Remove rows with any NAs
  complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  cat("Data for LSEI:", nrow(founder_matrix_clean), "complete SNPs\n")
  
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
    } else {
      cat("✗ LSEI failed with error code:", lsei_result$IsError, "\n")
      estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
    }
  }, error = function(e) {
    cat("✗ LSEI error:", e$message, "\n")
    estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
  })
} else {
  cat("Skipping LSEI - founders not distinguishable\n")
}

cat("\n=== FINAL RESULT ===\n")
cat("final_window_size:", final_window_size, "bp (", final_window_size/1000, "kb)\n")
cat("estimate_OK:", estimate_OK, "\n")
cat("n_snps:", final_n_snps, "\n")

if (estimate_OK == 1) {
  cat("✓ SUCCESS: Adaptive algorithm found optimal window size!\n")
  cat("✓ All", length(founders), "founders can be distinguished at", final_window_size/1000, "kb\n")
  cat("✓ LSEI estimation successful - haplotype frequencies available\n")
} else if (estimate_OK == 0) {
  cat("✗ FAILURE: Could not distinguish all", length(founders), "founders\n")
  cat("✗ Even at maximum window size of", max(window_sizes)/1000, "kb\n")
  cat("✗ Best result:", max(snps_tracking$n_groups), "groups\n")
} else {
  cat("✗ LSEI estimation failed - no haplotype frequencies available\n")
}

cat("\n=== ADAPTIVE WINDOW TEST COMPLETE ===\n")
