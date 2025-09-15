#!/usr/bin/env Rscript

# Test script to reproduce Friday night results for one position
# Usage: Rscript scripts/ErrMatrix/test_single_position.R

# Load required libraries
library(tidyverse)
library(limSolve)
library(MASS)
library(purrr)

# Read the Friday night version of BASE_VAR_WIDE.R and extract just the functions
base_var_wide_content <- readLines("scripts/ErrMatrix/BASE_VAR_WIDE.R")

# Find the function definitions (lines that start with function names)
function_lines <- grep("^[a-zA-Z_][a-zA-Z0-9_]*\\s*<-\\s*function", base_var_wide_content)

# Source the actual Friday night working code (BASE_VAR_WIDE.R)
# Override interactive() to avoid command line parsing
old_interactive <- interactive
interactive <- function() TRUE  # Override interactive() to return TRUE
source("scripts/ErrMatrix/BASE_VAR_WIDE.R")
interactive <- old_interactive  # Restore original

# Test parameters
chr <- "chr3R"
testing_position <- 19610000
sample_name <- "Rep01_W_F"
param_file <- "helpfiles/ZINC2_haplotype_parameters.R"
output_dir <- "process/ZINC2"

cat("=== TESTING SINGLE POSITION ===\n")
cat("Chromosome:", chr, "\n")
cat("Position:", testing_position, "\n")
cat("Sample:", sample_name, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load parameters
source(param_file, local = TRUE)
founders <- get("founders")
parameter <- get("h_cutoff")  # The parameter file defines h_cutoff, not parameter
method <- "adaptive"  # Set method directly since it's not in the parameter file

cat("Loaded parameters:\n")
cat("  Founders:", paste(founders, collapse=", "), "\n")
cat("  Parameter:", parameter, "\n")
cat("  Method:", method, "\n\n")

# Load RefAlt data using BASE_VAR_WIDE.R function
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
cat("Loading RefAlt data from:", refalt_file, "\n")

if (!file.exists(refalt_file)) {
  stop("RefAlt file not found: ", refalt_file)
}

# Use the process_refalt_data function from BASE_VAR_WIDE.R
df3 <- process_refalt_data(refalt_file, founders)
cat("Processed", nrow(df3), "rows for", ncol(df3)-1, "samples\n\n")

# Create subset for this position (500kb window)
window_size <- 500000
df4 <- df3 %>%
  dplyr::filter(POS >= testing_position - window_size/2 & 
                POS <= testing_position + window_size/2)

cat("Subsetted data for position", testing_position, ":\n")
cat("  Window size:", window_size, "bp\n")
cat("  SNPs in window:", nrow(df4), "\n")
cat("  Columns:", paste(names(df4), collapse=", "), "\n\n")

# Save the subsetted data
subset_file <- paste0("test_subset_", chr, "_", testing_position, "_", sample_name, ".RDS")
saveRDS(df4, subset_file)
cat("Saved subsetted data to:", subset_file, "\n\n")

# Run the estimator
cat("Running estimator...\n")
result <- estimate_haplotypes_list_format(
  pos = testing_position,
  sample_name = sample_name,
  df3 = df4,
  founders = founders,
  h_cutoff = parameter,
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

# Save the result
result_file <- paste0("test_result_", chr, "_", testing_position, "_", sample_name, ".RDS")
saveRDS(result, result_file)
cat("\nSaved result to:", result_file, "\n")

# Load Friday night production results for comparison
cat("\n=== LOADING FRIDAY NIGHT PRODUCTION RESULTS ===\n")
prod_file <- file.path(output_dir, "haplotype_results_list_format", "adapt_h4", paste0("R.haps.", chr, ".out.rds"))
cat("Loading production results from:", prod_file, "\n")

if (file.exists(prod_file)) {
  prod_data <- readRDS(prod_file)
  
  # Find the matching position and sample
  prod_entry <- prod_data %>%
    dplyr::filter(CHROM == chr, pos == testing_position) %>%
    tidyr::unnest(c(sample, Groups, Haps, Err, Names)) %>%
    dplyr::filter(sample == sample_name)
  
  if (nrow(prod_entry) > 0) {
    prod_groups <- prod_entry$Groups[[1]]
    prod_haps <- prod_entry$Haps[[1]]
    prod_err <- prod_entry$Err[[1]]
    prod_names <- prod_entry$Names[[1]]
    
    cat("\n=== FRIDAY NIGHT PRODUCTION RESULT ===\n")
    cat("Groups:", paste(prod_groups, collapse=","), "\n")
    cat("Haplotypes:\n")
    print(prod_haps)
    cat("Error matrix:\n")
    print(prod_err)
    cat("Names:", paste(prod_names, collapse=","), "\n")
    
    # Compare results
    cat("\n=== COMPARISON ===\n")
    groups_match <- identical(sort(result$Groups), sort(prod_groups))
    names_match <- identical(sort(result$Names), sort(prod_names))
    
    cat("Groups match:", groups_match, "\n")
    cat("Names match:", names_match, "\n")
    
    if (is.matrix(result$Err) && is.matrix(prod_err)) {
      # Align matrices by founder names
      debug_founders <- rownames(result$Err)
      prod_founders <- rownames(prod_err)
      
      if (all(sort(debug_founders) == sort(prod_founders))) {
        # Reorder debug matrix to match production order
        debug_err_aligned <- result$Err[prod_founders, prod_founders]
        
        err_diff <- abs(debug_err_aligned - prod_err)
        max_diff <- max(err_diff, na.rm = TRUE)
        sum_diff <- sum(err_diff, na.rm = TRUE)
        
        cat("Error matrix max difference:", signif(max_diff, 6), "\n")
        cat("Error matrix sum difference:", signif(sum_diff, 6), "\n")
        cat("Error matrices match:", max_diff < 1e-10, "\n")
        
        if (max_diff > 1e-10) {
          cat("\nDifference matrix (test - production):\n")
          print(signif(debug_err_aligned - prod_err, 4))
        }
      } else {
        cat("Founder names don't match - cannot compare error matrices\n")
        cat("Test founders:", paste(debug_founders, collapse=","), "\n")
        cat("Prod founders:", paste(prod_founders, collapse=","), "\n")
      }
    } else {
      cat("Cannot compare error matrices - one or both are not matrices\n")
    }
    
  } else {
    cat("Position/sample not found in Friday night production output\n")
  }
} else {
  cat("Friday night production file not found:", prod_file, "\n")
}

cat("\nDone.\n")
