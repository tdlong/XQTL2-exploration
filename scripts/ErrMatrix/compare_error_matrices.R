#!/usr/bin/env Rscript

# Script to compare error matrices between production and simple estimator test
# Usage: Rscript scripts/ErrMatrix/compare_error_matrices.R

library(tidyverse)

# Parameters
chr <- "chr3R"
position <- 19610000
sample_name <- "Rep01_W_F"
prod_file <- "process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds"
test_file <- "position_data_chr3R_19610000_Rep01_W_F.RDS"

cat("=== COMPARING ERROR MATRICES ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", position, "\n")
cat("Sample:", sample_name, "\n")
cat("Production file:", prod_file, "\n")
cat("Test file:", test_file, "\n\n")

# Load production results
cat("Loading production results...\n")
if (!file.exists(prod_file)) {
  stop("Production file not found: ", prod_file)
}
prod_data <- readRDS(prod_file)

# Extract production error matrix for this position and sample
prod_entry <- prod_data %>%
  dplyr::filter(CHROM == chr, pos == position) %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names)) %>%
  dplyr::filter(sample == sample_name)

if (nrow(prod_entry) == 0) {
  stop("Position/sample not found in production results")
}

prod_groups <- prod_entry$Groups[[1]]
prod_haps <- prod_entry$Haps[[1]]
prod_err <- prod_entry$Err[[1]]
prod_names <- prod_entry$Names[[1]]

cat("Production results:\n")
cat("  Groups:", paste(prod_groups, collapse=","), "\n")
cat("  Names:", paste(prod_names, collapse=","), "\n")
cat("  Error matrix dimensions:", nrow(prod_err), "x", ncol(prod_err), "\n\n")

# Load test results
cat("Loading test results...\n")
if (!file.exists(test_file)) {
  stop("Test file not found: ", test_file)
}

# Load BASE_VAR_WIDE.R functions
old_interactive <- interactive
interactive <- function() TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive

# Load test data and run estimator
payload <- readRDS(test_file)
df4 <- payload$df4
testing_position <- payload$testing_position
founders <- payload$founders
h_cutoff <- payload$h_cutoff
method <- payload$method
chr <- payload$chr

cat("Running test estimator...\n")
test_result <- est_haps_var(
  testing_position = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = h_cutoff,
  method = method,
  chr = chr,
  verbose = 2
)

test_groups <- test_result$Groups
test_haps <- test_result$Haps
test_err <- test_result$Err
test_names <- test_result$Names

cat("Test results:\n")
cat("  Groups:", paste(test_groups, collapse=","), "\n")
cat("  Names:", paste(test_names, collapse=","), "\n")
cat("  Error matrix dimensions:", nrow(test_err), "x", ncol(test_err), "\n\n")

# Compare results
cat("=== COMPARISON ===\n")
groups_match <- identical(sort(test_groups), sort(prod_groups))
names_match <- identical(sort(test_names), sort(prod_names))

cat("Groups match:", groups_match, "\n")
cat("Names match:", names_match, "\n")

if (is.matrix(test_err) && is.matrix(prod_err)) {
  # Align matrices by founder names
  test_founders <- rownames(test_err)
  prod_founders <- rownames(prod_err)
  
  if (all(sort(test_founders) == sort(prod_founders))) {
    # Reorder test matrix to match production order
    test_err_aligned <- test_err[prod_founders, prod_founders]
    
    err_diff <- abs(test_err_aligned - prod_err)
    max_diff <- max(err_diff, na.rm = TRUE)
    sum_diff <- sum(err_diff, na.rm = TRUE)
    
    cat("Error matrix max absolute difference:", signif(max_diff, 6), "\n")
    cat("Error matrix sum absolute difference:", signif(sum_diff, 6), "\n")
    cat("Error matrices match (tolerance 1e-10):", max_diff < 1e-10, "\n")
    
    if (max_diff > 1e-10) {
      cat("\nDifference matrix (test - production):\n")
      print(signif(test_err_aligned - prod_err, 4))
    }
  } else {
    cat("Founder names don't match - cannot compare error matrices\n")
    cat("Test founders:", paste(test_founders, collapse=","), "\n")
    cat("Prod founders:", paste(prod_founders, collapse=","), "\n")
  }
} else {
  cat("Cannot compare error matrices - one or both are not matrices\n")
}

cat("\nDone.\n")
