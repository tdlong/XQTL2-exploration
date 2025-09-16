#!/usr/bin/env Rscript

# Script to debug data consistency between production and test runs
# Usage: Rscript scripts/ErrMatrix/debug_data_consistency.R

library(tidyverse)

# Parameters
chr <- "chr3R"
position <- 19610000
sample_name <- "Rep01_W_F"

cat("=== DEBUGGING DATA CONSISTENCY ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", position, "\n")
cat("Sample:", sample_name, "\n\n")

# Load the test data
test_file <- "position_data_chr3R_19610000_Rep01_W_F.RDS"
payload <- readRDS(test_file)
df4 <- payload$df4
testing_position <- payload$testing_position
founders <- payload$founders
h_cutoff <- payload$h_cutoff
method <- payload$method
chr <- payload$chr

cat("Test data summary:\n")
cat("  Position:", testing_position, "\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Position range:", min(df4$POS), "to", max(df4$POS), "\n")
cat("  Founders:", paste(founders, collapse=", "), "\n")
cat("  h_cutoff:", h_cutoff, "\n")
cat("  Method:", method, "\n\n")

# Check if the sample exists in the data
if (sample_name %in% names(df4)) {
  cat("Sample", sample_name, "found in data\n")
  sample_data <- df4[[sample_name]]
  cat("  Sample data range:", min(sample_data, na.rm=TRUE), "to", max(sample_data, na.rm=TRUE), "\n")
  cat("  Sample data NAs:", sum(is.na(sample_data)), "\n")
} else {
  cat("Sample", sample_name, "NOT found in data\n")
  cat("Available samples:", paste(names(df4)[!names(df4) %in% c("POS", founders)], collapse=", "), "\n")
}

# Check founder data
cat("\nFounder data summary:\n")
for (founder in founders) {
  if (founder %in% names(df4)) {
    founder_data <- df4[[founder]]
    cat("  ", founder, ": range", min(founder_data, na.rm=TRUE), "to", max(founder_data, na.rm=TRUE), 
        ", NAs:", sum(is.na(founder_data)), "\n")
  } else {
    cat("  ", founder, ": NOT FOUND\n")
  }
}

# Check data around the target position
cat("\nData around target position", position, ":\n")
nearby_data <- df4 %>%
  filter(POS >= position - 1000 & POS <= position + 1000) %>%
  arrange(POS)

cat("  SNPs within 1kb of target:", nrow(nearby_data), "\n")
if (nrow(nearby_data) > 0) {
  cat("  Closest positions:", paste(nearby_data$POS[1:min(5, nrow(nearby_data))], collapse=", "), "\n")
}

cat("\nDone.\n")
