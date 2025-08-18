#!/usr/bin/env Rscript

# Test Fixed Window Algorithm
# Based on the proven adaptive window test, but single window size
# This tests the CORE algorithm that both production scripts should use

library(tidyverse)
library(limSolve)

cat("=== FIXED WINDOW ALGORITHM TEST ===\n\n")

# 1. Load parameters (same as adaptive test)
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and transform data (EXACT same as adaptive test)
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

# Filter for high-quality SNPs (EXACT same as adaptive test)
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("✓ Data ready:", nrow(df3), "rows\n")

# 3. Test the fixed window algorithm
cat("\n3. Testing fixed window algorithm...\n")

# Test parameters
test_pos <- 15000000
test_window_size <- 50000  # Fixed 50kb window
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test window size:", test_window_size, "bp\n")
cat("✓ Test sample:", test_sample, "\n\n")

# Create window around test position
window_start <- test_pos - test_window_size/2
window_end <- test_pos + test_window_size/2
cat("Window range:", window_start, "to", window_end, "bp\n")

# Get SNPs in window (long format first)
window_snps_long <- df3 %>%
  filter(CHROM == "chr2R" &
         POS > window_start &
         POS < window_end &
         (name %in% founders | name == test_sample))

# Convert to wide format (rows = positions, columns = sample + founders)
wide_data <- window_snps_long %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("SNPs in window:", nrow(wide_data), "positions\n")

# Quality filter: keep rows where ALL founders are fixed (freq < 0.03 OR freq > 0.97)
quality_filtered <- wide_data %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  )

cat("After quality filter:", nrow(quality_filtered), "positions\n")

# Check if we have enough data
if (nrow(quality_filtered) < 10) {
  cat("✗ Insufficient data for estimation\n")
  quit()
}

# Get founder matrix and sample frequencies from wide format
founder_matrix <- quality_filtered %>%
  select(all_of(founders)) %>%
  as.matrix()

# Get sample frequencies
sample_freqs <- quality_filtered %>%
  pull(!!test_sample)

cat("Sample frequencies available:", length(sample_freqs), "positions\n")

# Filter for non-NA values
valid_positions <- !is.na(sample_freqs)
sample_freqs <- sample_freqs[valid_positions]
founder_matrix <- founder_matrix[valid_positions, ]

cat("Final valid positions:", sum(valid_positions), "out of", length(valid_positions), "\n")

if (nrow(founder_matrix) < 10) {
  cat("✗ Insufficient valid data for estimation\n")
  quit()
}

# Run hierarchical clustering (SAME as adaptive test)
founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
distances <- dist(t(founder_matrix_clean), method = "euclidean")
hclust_result <- hclust(distances, method = "ward.D2")

# Use h_cutoff from parameter file
groups <- cutree(hclust_result, h = h_cutoff)
n_groups <- length(unique(groups))

cat("Founder groups:", n_groups, "\n")
cat("Group assignments:", paste(groups, collapse=", "), "\n")

# Run LSEI (SAME as adaptive test, but no constraint accumulation)
n_founders <- ncol(founder_matrix_clean)
E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
F <- 1.0

cat("Running LSEI with", n_founders, "founders...\n")

# Solve constrained least squares
tryCatch({
  result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs[complete.cases(founder_matrix)], 
                          E = E, F = F, 
                          G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
  
  if (result$IsError == 0) {
    cat("✓ LSEI successful!\n")
    cat("Estimated frequencies:", paste(sprintf("%.4f", result$X), collapse=", "), "\n")
    cat("Sum of frequencies:", sum(result$X), "\n")
    
    # Show which founders are in which groups
    for (group_id in unique(groups)) {
      group_founders <- founders[groups == group_id]
      group_freq <- sum(result$X[groups == group_id])
      cat("Group", group_id, ":", paste(group_founders, collapse="+"), "=", round(group_freq, 4), "\n")
    }
    
  } else {
    cat("✗ LSEI failed with error code:", result$IsError, "\n")
  }
}, error = function(e) {
  cat("✗ LSEI error:", e$message, "\n")
})

cat("\n=== FIXED WINDOW TEST COMPLETE ===\n")
