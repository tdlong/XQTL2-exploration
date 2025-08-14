#!/usr/bin/env Rscript

# Check Haplotype Data
# Examine the actual frequency values for this sample

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== Checking Haplotype Data ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n\n")

# Load haplotype results
h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
haplotype_results <- read_rds(haplotype_file) %>%
  filter(h_cutoff == !!h_cutoff)

cat("Total haplotype results:", nrow(haplotype_results), "\n")
cat("Samples:", paste(unique(haplotype_results$sample), collapse = ", "), "\n\n")

# Filter to sample and euchromatin
euchromatin_start <- 5398184
euchromatin_end <- 24684540

sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name, 
         pos >= euchromatin_start, 
         pos <= euchromatin_end)

cat("Sample haplotypes in euchromatin:", nrow(sample_haplotypes), "\n")

if (nrow(sample_haplotypes) == 0) {
  cat("❌ NO haplotypes for this sample in euchromatin!\n")
  quit(status = 1)
}

# Check frequency values
cat("\n=== Frequency Analysis ===\n")
freq_summary <- sample_haplotypes %>%
  summarise(
    total_rows = n(),
    na_freq = sum(is.na(freq)),
    non_na_freq = sum(!is.na(freq)),
    na_percent = round(na_freq / total_rows * 100, 1)
  )

print(freq_summary)

# Check a few specific positions
cat("\n=== Sample Positions ===\n")
test_positions <- c(5400000, 5410000, 5420000)

for (pos in test_positions) {
  pos_data <- sample_haplotypes %>% filter(pos == !!pos)
  cat("Position", pos, ":\n")
  if (nrow(pos_data) > 0) {
    for (i in 1:nrow(pos_data)) {
      cat("  Founder:", pos_data$founder[i], "Freq:", pos_data$freq[i], "\n")
    }
  } else {
    cat("  No data for this position\n")
  }
  cat("\n")
}

# Check all founders for a specific position
cat("=== All Founders at Position 5400000 ===\n")
pos_5400000 <- sample_haplotypes %>% filter(pos == 5400000)
if (nrow(pos_5400000) > 0) {
  for (i in 1:nrow(pos_5400000)) {
    cat("Founder:", pos_5400000$founder[i], "Freq:", pos_5400000$freq[i], "\n")
  }
} else {
  cat("No data for position 5400000\n")
}

# Check if this sample has any non-NA frequencies anywhere
cat("\n=== Overall Sample Check ===\n")
all_sample_data <- haplotype_results %>% filter(sample == sample_name)
cat("Total rows for sample:", nrow(all_sample_data), "\n")
cat("Non-NA frequencies:", sum(!is.na(all_sample_data$freq)), "\n")
cat("NA frequencies:", sum(is.na(all_sample_data$freq)), "\n")

if (sum(!is.na(all_sample_data$freq)) == 0) {
  cat("❌ ALL frequencies are NA for this sample!\n")
  cat("This explains why interpolation fails.\n")
} else {
  cat("✓ Some non-NA frequencies exist\n")
}

cat("\n=== Check Complete ===\n")
