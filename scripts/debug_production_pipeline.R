#!/usr/bin/env Rscript

# Debug Production Pipeline - Step by Step Data Inspection
# This script shows the data at each step to identify where the production pipeline fails

cat("=== DEBUG PRODUCTION PIPELINE STEP BY STEP ===\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load REFALT data (same as production script)
cat("\n2. Loading REFALT data...\n")
filein <- "process/JUICE/RefAlt.chr2R.txt"
df <- read.table(filein, header = TRUE, sep = "\t", fill = TRUE)
cat("✓ Raw data loaded:", nrow(df), "rows,", ncol(df), "columns\n")
cat("✓ POS column range:", range(df$POS, na.rm = TRUE), "\n")
cat("✓ POS column class:", class(df$POS), "\n")
cat("✓ First few POS values:", head(df$POS, 10), "\n")

# 3. Transform to frequencies (same as production script)
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
cat("✓ POS column range after transform:", range(df2$POS, na.rm = TRUE), "\n")
cat("✓ POS column class after transform:", class(df2$POS), "\n")

# 4. Filter for high-quality SNPs (same as production script)
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
cat("✓ Good SNPs POS range:", range(good_snps$POS, na.rm = TRUE), "\n")

# 5. Subset to high-quality SNPs (same as production script)
cat("\n5. Subsetting to high-quality SNPs...\n")
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

cat("✓ Subset data:", nrow(df3), "rows\n")
cat("✓ Subset POS range:", range(df3$POS, na.rm = TRUE), "\n")

# 6. Get non-founder samples (same as production script)
cat("\n6. Getting non-founder samples...\n")
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("✓ Non-founder samples:", length(non_founder_samples), "\n")
cat("✓ Sample names:", paste(non_founder_samples, collapse=", "), "\n")

# 7. Calculate chromosome length (same as production script)
cat("\n7. Calculating chromosome length...\n")
chromosome_length <- max(df$POS, na.rm = TRUE)
cat("✓ Raw chromosome length:", chromosome_length, "\n")
cat("✓ Is finite:", is.finite(chromosome_length), "\n")
cat("✓ Is NA:", is.na(chromosome_length), "\n")

# 8. Define scan positions (same as production script)
cat("\n8. Defining scan positions...\n")
scan_start <- 500000
scan_end <- chromosome_length - 500000
cat("✓ Scan range:", scan_start, "to", scan_end, "\n")
cat("✓ Scan end is finite:", is.finite(scan_end), "\n")

if (is.finite(scan_end) && scan_end > scan_start) {
  scan_positions <- seq(scan_start, scan_end, by = 10000)
  cat("✓ Scan positions created:", length(scan_positions), "positions\n")
  cat("✓ First few scan positions:", head(scan_positions, 5), "\n")
} else {
  cat("✗ Cannot create scan positions - invalid range\n")
}

cat("\n=== PIPELINE DEBUG COMPLETE ===\n")
