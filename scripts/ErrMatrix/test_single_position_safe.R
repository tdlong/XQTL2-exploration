l runs t#!/usr/bin/env Rscript

# Test script to reproduce Friday night results for ONE position only
# This runs BASE_VAR_WIDE.R but only for position 19610000
# Saves to a different directory to avoid overwriting Friday night files

# Set up the same arguments as Friday night but with different output directory
chr <- "chr3R"
method <- "adaptive"
parameter <- 4
output_dir <- "process/ZINC2_test"  # Different directory to avoid overwriting
param_file <- "helpfiles/ZINC2_haplotype_parameters.R"
target_position <- 19610000

cat("=== TESTING SINGLE POSITION REPRODUCTION ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Target position:", target_position, "\n\n")

# Create test output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "haplotype_results_list_format", "adapt_h4"), recursive = TRUE, showWarnings = FALSE)

# Copy RefAlt file to test directory
refalt_source <- "process/ZINC2/RefAlt.chr3R.txt"
refalt_dest <- file.path(output_dir, "RefAlt.chr3R.txt")
if (file.exists(refalt_source)) {
  file.copy(refalt_source, refalt_dest, overwrite = TRUE)
  cat("Copied RefAlt file to test directory\n")
} else {
  stop("RefAlt file not found: ", refalt_source)
}

# Modify BASE_VAR_WIDE.R to only process our target position
# We'll do this by creating a modified version that limits positions
cat("Creating modified BASE_VAR_WIDE.R for single position test...\n")

# Read the original BASE_VAR_WIDE.R
base_content <- readLines("scripts/ErrMatrix/BASE_VAR_WIDE.R")

# Find the line that defines all_positions and replace it
# Look for: all_positions <- seq(scan_start, scan_end, by = step)
pos_line <- grep("all_positions.*seq.*scan_start.*scan_end", base_content)
if (length(pos_line) > 0) {
  base_content[pos_line] <- paste0("all_positions <- c(", target_position, ")  # Only process target position")
  cat("Modified all_positions to only include position", target_position, "\n")
} else {
  cat("Warning: Could not find all_positions line to modify\n")
}

# Write modified version
modified_file <- "scripts/ErrMatrix/BASE_VAR_WIDE_single.R"
writeLines(base_content, modified_file)

# Set up command line arguments for the modified script
args <- c(chr, method, parameter, output_dir, param_file, "--nonverbose")

# Override commandArgs to simulate the Friday night call
commandArgs <- function(trailingOnly = TRUE) args

# Source and run the modified BASE_VAR_WIDE.R
cat("Running modified BASE_VAR_WIDE.R for single position...\n")
source(modified_file)

cat("\n=== COMPARISON WITH FRIDAY NIGHT RESULTS ===\n")

# Load Friday night production results
prod_file <- "process/ZINC2/haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds"
if (file.exists(prod_file)) {
  cat("Loading Friday night production results...\n")
  prod_data <- readRDS(prod_file)
  
  # Find the matching position and sample
  prod_entry <- prod_data %>%
    dplyr::filter(CHROM == chr, pos == target_position) %>%
    tidyr::unnest(c(sample, Groups, Haps, Err, Names)) %>%
    dplyr::filter(sample == "Rep01_W_F")
  
  if (nrow(prod_entry) > 0) {
    prod_groups <- prod_entry$Groups[[1]]
    prod_haps <- prod_entry$Haps[[1]]
    prod_err <- prod_entry$Err[[1]]
    prod_names <- prod_entry$Names[[1]]
    
    cat("Friday night production result for position", target_position, ":\n")
    cat("Groups:", paste(prod_groups, collapse=","), "\n")
    cat("Haplotypes:\n")
    print(prod_haps)
    cat("Error matrix:\n")
    print(prod_err)
    cat("Names:", paste(prod_names, collapse=","), "\n")
  } else {
    cat("Position", target_position, "not found in Friday night production results\n")
  }
} else {
  cat("Friday night production file not found:", prod_file, "\n")
}

# Load our test results
test_file <- file.path(output_dir, "haplotype_results_list_format", "adapt_h4", "R.haps.chr3R.out.rds")
if (file.exists(test_file)) {
  cat("\nLoading test results...\n")
  test_data <- readRDS(test_file)
  
  # Find the matching position and sample
  test_entry <- test_data %>%
    dplyr::filter(CHROM == chr, pos == target_position) %>%
    tidyr::unnest(c(sample, Groups, Haps, Err, Names)) %>%
    dplyr::filter(sample == "Rep01_W_F")
  
  if (nrow(test_entry) > 0) {
    test_groups <- test_entry$Groups[[1]]
    test_haps <- test_entry$Haps[[1]]
    test_err <- test_entry$Err[[1]]
    test_names <- test_entry$Names[[1]]
    
    cat("Test result for position", target_position, ":\n")
    cat("Groups:", paste(test_groups, collapse=","), "\n")
    cat("Haplotypes:\n")
    print(test_haps)
    cat("Error matrix:\n")
    print(test_err)
    cat("Names:", paste(test_names, collapse=","), "\n")
    
    # Compare results
    if (nrow(prod_entry) > 0) {
      cat("\n=== COMPARISON ===\n")
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
          
          cat("Error matrix max difference:", signif(max_diff, 6), "\n")
          cat("Error matrix sum difference:", signif(sum_diff, 6), "\n")
          cat("Error matrices match:", max_diff < 1e-10, "\n")
        } else {
          cat("Founder names don't match - cannot compare error matrices\n")
        }
      }
    }
  } else {
    cat("Position", target_position, "not found in test results\n")
  }
} else {
  cat("Test results file not found:", test_file, "\n")
}

# Clean up
unlink(modified_file)
cat("\nDone.\n")
