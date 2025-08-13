#!/usr/bin/env Rscript

# Quick Haplotype Debug Script
# Fast diagnostic tool to examine haplotype data without running full pipeline

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Set parameters
chr <- "chr2R"
output_dir <- "process/JUICE"
estimator <- "adaptive_h6"
sample_name <- "GJ_3_1"

cat("=== Quick Haplotype Debug ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample:", sample_name, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load haplotype results
cat("Loading haplotype results...\n")
if (grepl("^adaptive_h", estimator)) {
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, paste0("adaptive_window_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file) %>%
    filter(h_cutoff == !!h_cutoff)
  cat("✓ Adaptive window results loaded for h_cutoff =", h_cutoff, "\n")
} else {
  stop("Only adaptive estimators supported for now")
}

cat("Total haplotype estimates:", nrow(haplotype_results), "\n")
cat("Samples:", paste(unique(haplotype_results$sample), collapse = ", "), "\n\n")

# Filter to specific sample
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

cat("=== Sample Data Analysis ===\n")
cat("Sample:", sample_name, "\n")
cat("Total haplotype estimates:", nrow(sample_haplotypes), "\n")

if (nrow(sample_haplotypes) == 0) {
  cat("❌ NO HAPLOTYPE ESTIMATES for this sample!\n")
  cat("This explains the failure.\n")
  quit(status = 1)
}

# Check position ranges
cat("Position range:", min(sample_haplotypes$pos), "-", max(sample_haplotypes$pos), "bp\n")

# Check euchromatin region
euchromatin_start <- 5398184
euchromatin_end <- 24684540

euchromatic_haplotypes <- sample_haplotypes %>%
  filter(pos >= euchromatin_start & pos <= euchromatin_end)

cat("Euchromatin region haplotypes:", nrow(euchromatic_haplotypes), "\n")

if (nrow(euchromatic_haplotypes) == 0) {
  cat("❌ NO HAPLOTYPES in euchromatin region!\n")
  cat("All haplotypes are outside the target region.\n")
  quit(status = 1)
}

# Check for NA frequencies
na_check <- euchromatic_haplotypes %>%
  filter(!is.na(freq)) %>%
  nrow()

cat("Non-NA frequencies:", na_check, "/", nrow(euchromatic_haplotypes), "\n")

if (na_check == 0) {
  cat("❌ ALL frequencies are NA!\n")
  cat("This explains the failure.\n")
  quit(status = 1)
}

# Check founder coverage
founders <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
founder_coverage <- euchromatic_haplotypes %>%
  filter(!is.na(freq)) %>%
  group_by(founder) %>%
  summarise(
    n_positions = n(),
    mean_freq = mean(freq, na.rm = TRUE),
    min_freq = min(freq, na.rm = TRUE),
    max_freq = max(freq, na.rm = TRUE)
  )

cat("\n=== Founder Coverage ===\n")
print(founder_coverage)

# Check position distribution
cat("\n=== Position Distribution ===\n")
position_summary <- euchromatic_haplotypes %>%
  filter(!is.na(freq)) %>%
  summarise(
    n_positions = n(),
    min_pos = min(pos),
    max_pos = max(pos),
    mean_spacing = (max(pos) - min(pos)) / (n() - 1)
  )

print(position_summary)

# Quick SNP overlap check
cat("\n=== SNP Overlap Check ===\n")
# Load a small sample of SNPs to check overlap
refalt_file <- file.path(output_dir, paste0("df3.", chr, ".RDS"))
if (file.exists(refalt_file)) {
  df2 <- read_rds(refalt_file)
  snp_positions <- df2 %>% distinct(POS) %>% pull(POS) %>% sort()
  
  # Check overlap with haplotype positions
  haplotype_positions <- euchromatic_haplotypes %>% 
    filter(!is.na(freq)) %>% 
    pull(pos) %>% 
    sort()
  
  overlap_count <- sum(snp_positions %in% haplotype_positions)
  cat("SNP positions:", length(snp_positions), "\n")
  cat("Haplotype positions:", length(haplotype_positions), "\n")
  cat("Overlapping positions:", overlap_count, "\n")
  cat("Overlap percentage:", round(overlap_count / length(snp_positions) * 100, 2), "%\n")
  
  if (overlap_count == 0) {
    cat("❌ NO OVERLAP between SNPs and haplotypes!\n")
    cat("This would cause interpolation to fail.\n")
  }
}

cat("\n=== Debug Complete ===\n")
