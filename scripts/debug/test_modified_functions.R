#!/usr/bin/env Rscript

# Test the modified functions that capture groups and error matrices
# This will call the modified working function and show the captured data

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/debug/test_modified_functions.R <chr> <param_file> <output_dir> <position>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
test_position <- as.numeric(args[4])

cat("=== TESTING MODIFIED FUNCTIONS WITH GROUPS AND ERROR MATRICES ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Source the MODIFIED functions (not the original ones)
source("scripts/debug/haplotype_estimation_functions_with_groups.R")

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

# Load REFALT data (same as existing working code)
refalt_data <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(refalt_data), "rows\n")

# Transform to frequencies (same as existing working code)
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

# Filter for high-quality SNPs (same as existing working code)
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

# Get valid SNPs for euchromatin (same as existing working code)
observed_data <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(refalt_processed, multiple = "all")

# Filter to euchromatin and only the samples we actually processed (same as existing working code)
observed_euchromatic <- observed_data %>%
  filter(POS >= euchromatin_start, POS <= euchromatin_end) %>%
  filter(name %in% names_in_bam)

cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Test the modified function on the first sample
test_sample <- names_in_bam[1]  # Use first sample
cat("=== TESTING MODIFIED FUNCTION ===\n")
cat("Test position:", test_position, "\n")
cat("Test sample:", test_sample, "\n\n")

# Call the modified function
modified_result <- estimate_haplotypes(test_position, test_sample, observed_euchromatic, 
                                     founders, h_cutoff, chr, verbose = 2)

cat("\n=== MODIFIED FUNCTION RESULTS ===\n")
cat("Estimate OK:", modified_result$estimate_OK, "\n")
cat("Final window size:", modified_result$final_window_size, "\n")
cat("Number of SNPs:", modified_result$n_snps, "\n")
cat("Haplotype frequencies:\n")
for (founder in founders) {
  cat(sprintf("  %s: %.3f\n", founder, modified_result[[founder]]))
}

# Show the captured groups and error matrix
cat("\n=== CAPTURED GROUPS AND ERROR MATRIX ===\n")
if (!is.null(modified_result$groups)) {
  cat("Groups:\n")
  print(modified_result$groups)
} else {
  cat("Groups: NULL\n")
}

if (!is.null(modified_result$error_matrix)) {
  cat("Error matrix dimensions:", dim(modified_result$error_matrix), "\n")
  if (!all(is.na(modified_result$error_matrix))) {
    cat("Error matrix (first few values):\n")
    print(modified_result$error_matrix[1:3, 1:3])
  } else {
    cat("Error matrix: All NAs\n")
  }
} else {
  cat("Error matrix: NULL\n")
}

cat("\n=== TEST COMPLETE ===\n")




