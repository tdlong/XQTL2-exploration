#!/usr/bin/env Rscript

# Simple script: load df4 data and call est_haps_var directly
# Usage: Rscript scripts/ErrMatrix/simple_estimator_test.R <data_file> <metadata_file> <sample_name>

library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Print simple df4 summary
print_df4_summary <- function(df4, founders, sample_name) {
  cat("--- df4 summary ---\n")
  cat("rows:", nrow(df4), "\n")
  cat("POS min/max/mean:", min(df4$POS, na.rm=TRUE), max(df4$POS, na.rm=TRUE), round(mean(df4$POS, na.rm=TRUE), 2), "\n")
  cols <- c(founders, sample_name)
  cols <- cols[cols %in% names(df4)]
  stats <- lapply(cols, function(cn) {
    v <- df4[[cn]]
    c(mean=round(mean(v, na.rm=TRUE),6), sd=round(sd(v, na.rm=TRUE),6), na=sum(is.na(v)))
  })
  stats_df <- as.data.frame(do.call(rbind, stats))
  rownames(stats_df) <- cols
  print(stats_df)
  cat("-------------------\n\n")
}

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

# Input summary for df4
print_df4_summary(df4, founders, sample_name)

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
