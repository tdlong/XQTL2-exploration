#!/usr/bin/env Rscript

# =============================================================================
# Debug Script for lsei Failing Case
# =============================================================================
# Tests the exact case that's producing invalid haplotype frequencies:
# Position: 510000, Sample: GJ_3_1, Window: 20kb
# Usage: Rscript scripts/debug_lsei_case.R

library(tidyverse)
library(limSolve)

# Set up the test case
test_chr <- "chr2R"
test_pos <- 510000
test_sample <- "GJ_3_1"
test_window_size <- 20000
test_mydir <- "process/JUICE"
test_parfile <- "helpfiles/JUICE/JUICE_haplotype_parameters.R"

cat("=== Debug lsei Failing Case ===\n")
cat("Chromosome:", test_chr, "\n")
cat("Position:", test_pos, "\n")
cat("Sample:", test_sample, "\n")
cat("Window size:", test_window_size, "bp\n")
cat("Parameter file:", test_parfile, "\n")
cat("Data directory:", test_mydir, "\n\n")

# Source the parameter file
if (file.exists(test_parfile)) {
  source(test_parfile)
  cat("✓ Parameter file loaded\n")
  cat("Founders:", paste(founders, collapse = ", "), "\n")
  cat("h_cutoff:", h_cutoff, "\n\n")
} else {
  stop("Parameter file not found: ", test_parfile)
}

# Load the data
cat("Loading data...\n")
refalt_file <- paste0(test_mydir, "/RefAlt.", test_chr, ".txt")
if (!file.exists(refalt_file)) {
  stop("REFALT file not found: ", refalt_file)
}

df <- read.table(refalt_file, header = TRUE)
cat("✓ REFALT data loaded:", nrow(df), "rows,", ncol(df), "columns\n")
cat("Columns:", paste(names(df), collapse = ", "), "\n\n")

# Transform REF/ALT counts to frequencies
cat("Transforming data...\n")
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

cat("✓ Data transformed\n")
cat("Unique names:", paste(unique(df2$name), collapse = ", "), "\n\n")

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

cat("✓ High-quality SNPs identified:", nrow(good_snps), "\n\n")

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Create window around test position
window_start <- max(0, test_pos - test_window_size/2)
window_end <- test_pos + test_window_size/2

cat("=== Window Analysis ===\n")
cat("Window:", window_start, "-", window_end, "bp\n")
cat("Window size:", test_window_size, "bp\n")

# Filter SNPs in this window
window_snps <- df3 %>%
  filter(CHROM == test_chr &
         POS > window_start &
         POS < window_end &
         (name %in% founders | name == test_sample)) %>%
  select(-c(CHROM, N)) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("✓ Window SNPs processed\n")
cat("SNPs in window:", nrow(window_snps), "\n")
cat("Columns after pivot_wider:", paste(names(window_snps), collapse = ", "), "\n\n")

# Check if we have the test sample
if (!(test_sample %in% names(window_snps))) {
  stop("Test sample", test_sample, "not found in window data")
}

# Get founder matrix and sample frequencies
founder_matrix <- window_snps %>% select(all_of(founders))
sample_freqs <- window_snps[[test_sample]]

cat("=== Matrix Analysis ===\n")
cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
cat("Sample frequencies length:", length(sample_freqs), "\n")
cat("Founder matrix columns:", paste(names(founder_matrix), collapse = ", "), "\n\n")

# Filter for non-NA values
valid_positions <- !is.na(sample_freqs)
sample_freqs <- sample_freqs[valid_positions]
founder_matrix <- founder_matrix[valid_positions, ]

cat("After filtering NAs:\n")
cat("Valid positions:", sum(valid_positions), "/", length(valid_positions), "\n")
cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
cat("Sample frequencies:", paste(sprintf("%.4f", sample_freqs), collapse = ", "), "\n\n")

if (nrow(founder_matrix) > 0) {
  # Convert to matrix for clustering
  founder_matrix <- as.matrix(founder_matrix)
  
  cat("=== Founder Matrix Details ===\n")
  cat("Matrix class:", class(founder_matrix), "\n")
  cat("Matrix dimensions:", dim(founder_matrix), "\n")
  cat("Matrix summary:\n")
  print(summary(founder_matrix))
  cat("\n")
  
  # Cluster founders based on similarity
  cat("=== Clustering Analysis ===\n")
  founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
  n_groups <- length(unique(founder_clusters))
  cat("Founder groups:", n_groups, "\n")
  for (group_id in unique(founder_clusters)) {
    group_founders <- names(founder_clusters[founder_clusters == group_id])
    cat("  Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
  }
  cat("\n")
  
  # Build constraint matrix
  cat("=== Constraint Matrix Construction ===\n")
  n_founders <- ncol(founder_matrix)
  E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
  F <- 1.0
  
  cat("Constraint matrix E:\n")
  print(E)
  cat("Constraint values F:\n")
  print(F)
  cat("\n")
  
  # Define bounds
  G <- diag(n_founders)  # Lower bounds (>= 0.0003)
  H <- matrix(rep(0.0003, n_founders))  # Minimum frequency
  
  cat("Bounds matrix G:\n")
  print(G)
  cat("Bounds vector H:\n")
  print(H)
  cat("\n")
  
  # Test lsei
  cat("=== Testing lsei ===\n")
  cat("Input matrices:\n")
  cat("A (founder_matrix) dimensions:", dim(founder_matrix), "\n")
  cat("B (sample_freqs) length:", length(sample_freqs), "\n")
  cat("E (constraints) dimensions:", dim(E), "\n")
  cat("F (constraint values) length:", length(F), "\n")
  cat("G (bounds) dimensions:", dim(G), "\n")
  cat("H (bound values) length:", length(H), "\n\n")
  
  tryCatch({
    result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                  G = G, H = H)
    
    cat("✓ lsei succeeded!\n")
    cat("Result X (haplotype frequencies):\n")
    for (i in seq_along(founders)) {
      cat("  ", founders[i], ":", sprintf("%.6f", result$X[i]), "\n")
    }
    cat("\n")
    
    cat("=== Validation ===\n")
    cat("Sum of frequencies:", sum(result$X), "\n")
    cat("Min frequency:", min(result$X), "\n")
    cat("Max frequency:", max(result$X), "\n")
    cat("All >= 0.0003:", all(result$X >= 0.0003), "\n")
    cat("Sum ≈ 1.0:", abs(sum(result$X) - 1.0) < 0.001, "\n")
    
  }, error = function(e) {
    cat("✗ lsei failed:", e$message, "\n")
    cat("Error details:", e, "\n")
  })
  
} else {
  cat("✗ No valid data for analysis\n")
}

cat("\n=== Debug Complete ===\n")
