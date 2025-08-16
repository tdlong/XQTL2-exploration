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
}

cat("\n=== TEST COMPLETE ===\n")
