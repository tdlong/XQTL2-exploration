#!/usr/bin/env Rscript

# Simple script: load df4 data and call est_haps_var directly
# Usage: Rscript scripts/ErrMatrix/simple_estimator_test.R <data_file> <metadata_file> <sample_name>

library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript scripts/ErrMatrix/simple_estimator_test.R <data_file> <sample_name>\n")
  cat("Example: Rscript scripts/ErrMatrix/simple_estimator_test.R position_data_chr3R_19610000_Rep01_W_F.RDS Rep01_W_F\n")
  quit(status = 1)
}

data_file <- args[1]
sample_name <- args[2]

cat("=== SIMPLE ESTIMATOR TEST ===\n")
cat("Data file:", data_file, "\n")
cat("Sample:", sample_name, "\n\n")

# Load the data
payload <- readRDS(data_file)
df4 <- payload$df4
testing_position <- payload$testing_position
founders <- payload$founders
h_cutoff <- payload$h_cutoff
method <- payload$method
chr <- payload$chr

cat("Loaded data:\n")
cat("  SNPs:", nrow(df4), "\n")
cat("  Founders:", paste(founders, collapse=", "), "\n")
cat("  Position:", testing_position, "\n")
cat("  h_cutoff:", h_cutoff, "\n\n")

# Load BASE_VAR_WIDE.R functions
cat("Loading BASE_VAR_WIDE.R functions...\n")
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Call est_haps_var directly
cat("Calling est_haps_var...\n")
result <- est_haps_var(
  testing_position = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = h_cutoff,
  method = method,
  chr = chr,
  verbose = 2
)

cat("\n=== RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse=","), "\n")
cat("Haplotypes:\n")
print(result$Haps)
cat("Error matrix:\n")
print(result$Err)
cat("Names:", paste(result$Names, collapse=","), "\n")

cat("\nDone.\n")
