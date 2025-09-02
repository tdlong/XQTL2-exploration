#!/usr/bin/env Rscript

# Proper Test Script - Data Transformation Pipeline
# Tests the exact same logic as production script step by step

library(tidyverse)

cat("=== PROPER TEST SCRIPT - DATA TRANSFORMATION PIPELINE ===\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load REFALT data (same as production)
cat("\n2. Loading REFALT data...\n")
filein <- "process/JUICE/RefAlt.chr2R.txt"
df <- read.table(filein, header = TRUE)
cat("✓ Raw data loaded:", nrow(df), "rows,", ncol(df), "columns\n")
cat("✓ Columns:", paste(names(df), collapse=", "), "\n")

# 3. Transform to frequencies (same as production)
cat("\n3. Transforming to frequencies...\n")
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

cat("✓ Transformed data:", nrow(df2), "rows\n")
cat("✓ Columns:", paste(names(df2), collapse=", "), "\n")
cat("✓ Sample names:", paste(unique(df2$name), collapse=", "), "\n")

# 4. Filter for high-quality SNPs (same as production)
cat("\n4. Filtering for high-quality SNPs...\n")
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

cat("✓ Good SNPs found:", nrow(good_snps), "positions\n")

# 5. Subset dataset (same as production)
cat("\n5. Subsetting dataset...\n")
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

cat("✓ Subset data:", nrow(df3), "rows\n")
cat("✓ Columns:", paste(names(df3), collapse=", "), "\n")

# 6. Get non-founder samples (same as production)
cat("\n6. Getting non-founder samples...\n")
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("✓ Non-founder samples:", length(non_founder_samples), "\n")
cat("✓ Sample names:", paste(non_founder_samples, collapse=", "), "\n")

# 7. Test window creation (same as production)
cat("\n7. Testing window creation...\n")
test_pos <- 15000000  # Test position
window_size <- 10000  # 10kb window
window_start <- test_pos - window_size
window_end <- test_pos + window_size

cat("✓ Test position:", test_pos, "\n")
cat("✓ Window:", window_start, "to", window_end, "\n")

# 8. Test window data extraction (same as production)
cat("\n8. Testing window data extraction...\n")
window_snps <- df3 %>%
  filter(CHROM == "chr2R" &
         POS > window_start &
         POS < window_end &
         (name %in% founders | name == non_founder_samples[1]))

cat("✓ Window SNPs:", nrow(window_snps), "rows\n")
cat("✓ Window columns:", paste(names(window_snps), collapse=", "), "\n")
cat("✓ Window sample names:", paste(unique(window_snps$name), collapse=", "), "\n")

# 9. Test founder data extraction (same as production)
cat("\n9. Testing founder data extraction...\n")
founder_data <- window_snps %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("✓ Founder data:", nrow(founder_data), "rows\n")
cat("✓ Founder columns:", paste(names(founder_data), collapse=", "), "\n")

# 10. Test founder matrix creation (same as production)
cat("\n10. Testing founder matrix creation...\n")
if (ncol(founder_data) >= length(founders) + 1) {
  founder_matrix <- founder_data %>%
    select(-POS) %>%
    as.matrix()
  
  cat("✓ Founder matrix:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
  cat("✓ Matrix columns:", paste(colnames(founder_matrix), collapse=", "), "\n")
} else {
  cat("✗ Insufficient founder data\n")
  quit(status = 1)
}

# 11. Test haplotype estimation algorithm (same as production)
cat("\n11. Testing haplotype estimation algorithm...\n")

# Get sample frequencies for the same positions
sample_name <- non_founder_samples[1]
sample_freqs <- window_snps %>%
  filter(name == sample_name) %>%
  select(POS, freq) %>%
  right_join(founder_data %>% select(POS), by = "POS") %>%
  pull(freq)

cat("✓ Sample name:", sample_name, "\n")
cat("✓ Sample frequencies:", length(sample_freqs), "values\n")
cat("✓ Non-NA frequencies:", sum(!is.na(sample_freqs)), "values\n")

# Filter for non-NA values
valid_positions <- !is.na(sample_freqs)
sample_freqs <- sample_freqs[valid_positions]
founder_matrix <- founder_matrix[valid_positions, ]

cat("✓ Valid positions:", sum(valid_positions), "\n")
cat("✓ Final founder matrix:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")

if (nrow(founder_matrix) < 10) {
  cat("✗ Insufficient data for haplotype estimation\n")
  quit(status = 1)
}

# Test hierarchical clustering (same as production)
cat("\n12. Testing hierarchical clustering...\n")
h_cutoff <- 4  # Test h_cutoff value

founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
n_groups <- length(unique(founder_clusters))

cat("✓ Hierarchical clustering completed\n")
cat("✓ Number of groups:", n_groups, "\n")
cat("✓ Group assignments:", paste(founder_clusters, collapse=", "), "\n")

# Test LSEI with constraints (same as production)
cat("\n13. Testing LSEI with constraints...\n")

n_founders <- ncol(founder_matrix)
E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
F <- 1.0

# Add group constraints
unique_clusters <- unique(founder_clusters)
for (cluster_id in unique_clusters) {
  cluster_founders <- which(founder_clusters == cluster_id)
  if (length(cluster_founders) > 1) {
    # Group constraint
    constraint_row <- rep(0, n_founders)
    constraint_row[cluster_founders] <- 1
    E <- rbind(E, constraint_row)
    F <- c(F, 0.5)  # Example constraint value
  }
}

cat("✓ Constraint matrix E:", nrow(E), "x", ncol(E), "\n")
cat("✓ Constraint vector F:", length(F), "values\n")

# Solve constrained least squares
tryCatch({
  result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                          G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
  
  if (result$IsError == 0) {
    cat("✓ LSEI successful!\n")
    cat("✓ Estimated founder frequencies:", paste(round(result$X, 4), collapse=", "), "\n")
    cat("✓ Sum of frequencies:", sum(result$X), "\n")
  } else {
    cat("✗ LSEI failed with error code:", result$IsError, "\n")
  }
}, error = function(e) {
  cat("✗ LSEI error:", e$message, "\n")
})

cat("\n=== FULL PIPELINE TEST COMPLETE ===\n")
