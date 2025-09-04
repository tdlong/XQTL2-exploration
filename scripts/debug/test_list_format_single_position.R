#!/usr/bin/env Rscript

# Test the actual list format haplotype estimation code on a single position
# This tests the estimate_haplotypes_list_format function from run_haplotype_estimation_list_format.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/debug/test_list_format_single_position.R <chr> <param_file> <output_dir> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])

cat("=== TESTING LIST FORMAT HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Define euchromatin boundaries (same as the main script)
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region:", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load observed data from REFALT files (same as the main script)
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load REFALT data (same as the main script)
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as the main script)
cat("Converting counts to frequencies...\n")
refalt_processed <- refalt_data %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(
    RefAlt = str_sub(lab, 1, 3),
    name = str_sub(lab, 5)
  ) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

# Filter for high-quality SNPs (same as the main script)
cat("Filtering for high-quality SNPs...\n")
good_snps <- refalt_processed %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(total_count == 0),
    not_fixed = sum(total_count != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Get valid SNPs for euchromatin (same as the main script)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed (same as the main script)
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Source the actual functions from the main script
cat("Loading the actual list format functions...\n")
source("scripts/production/run_haplotype_estimation_list_format.R")

# Test on first sample
test_sample <- names_in_bam[1]
cat("Testing on sample:", test_sample, "\n\n")

# Test the actual estimate_haplotypes_list_format function
cat("=== TESTING estimate_haplotypes_list_format FUNCTION ===\n")
result <- estimate_haplotypes_list_format(test_position, test_sample, observed_euchromatic, 
                                        founders, h_cutoff, "adaptive", 
                                        NULL, chr, verbose = 1)

cat("\n=== RESULTS ===\n")
cat("Result structure:\n")
str(result, max.level = 2)

cat("\nDetailed results:\n")
cat("CHROM:", result$CHROM, "\n")
cat("Position:", result$pos, "\n")
cat("Sample:", result$sample[[1]], "\n")

# Check groups
groups <- result$Groups[[1]][[1]]
cat("Groups:", paste(groups, collapse = ", "), "\n")
cat("Number of unique groups:", length(unique(groups)), "\n")
cat("Is 1:8?", all(sort(unique(groups)) == 1:8), "\n")

# Check haplotypes
haps <- result$Haps[[1]][[1]]
cat("Haplotypes:\n")
for (i in 1:length(haps)) {
  cat("  ", founders[i], ":", round(haps[i], 4), "\n")
}
cat("Sum:", round(sum(haps, na.rm = TRUE), 4), "\n")
cat("All NA?", all(is.na(haps)), "\n")

# Check error matrix
err <- result$Err[[1]][[1]]
cat("Error matrix dimensions:", nrow(err), "x", ncol(err), "\n")
cat("Error matrix all NA?", all(is.na(err)), "\n")
if (!all(is.na(err))) {
  cat("Error matrix has values:", sum(!is.na(err)), "/", length(err), "\n")
}

cat("\n=== TEST COMPLETE ===\n")
