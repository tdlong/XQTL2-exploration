#!/usr/bin/env Rscript

# Run Haplotype Estimation with List Format Output
# Focuses on adaptive_h4 and smooth_h4 methods
# Saves results to new hap_list_results directory

library(tidyverse)
library(limSolve)

# Source the new estimator function
source("scripts/production/haplotype_estimation_list_format.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript scripts/production/run_haplotype_estimation_list_format.R <chr> <param_file> <output_dir> <method>")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
method <- args[4]

# Validate method
if (!method %in% c("adaptive_h4", "smooth_h4")) {
  stop("Method must be 'adaptive_h4' or 'smooth_h4'")
}

cat("=== HAPLOTYPE ESTIMATION WITH LIST FORMAT ===\n")
cat("Chromosome:", chr, "\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Method:", method, "\n\n")

# Load parameters
source(param_file)

# Load observed data
observed_file <- file.path(output_dir, "observed_euchromatic", paste0("observed_euchromatic_", chr, ".RDS"))
if (!file.exists(observed_file)) {
  stop("Observed data file not found: ", observed_file)
}

observed_euchromatic <- readRDS(observed_file)
cat("✓ Observed data loaded:", nrow(observed_euchromatic), "rows\n")
cat("Samples:", paste(unique(observed_euchromatic$name), collapse = ", "), "\n\n")

# Create new results directory
hap_list_results_dir <- file.path(output_dir, "hap_list_results")
dir.create(hap_list_results_dir, showWarnings = FALSE, recursive = TRUE)

# Define positions to process
positions <- seq(euchromatin_start, euchromatin_end, by = 10000)
samples <- names_in_bam

cat("Processing", length(positions), "positions for", length(samples), "samples\n")
cat("Method:", method, "\n\n")

if (method == "adaptive_h4") {
  # Run adaptive_h4 with list format
  cat("Running adaptive_h4 estimation...\n")
  
  list_format_data <- positions %>%
    map_dfr(function(pos) {
      if (pos %% 100000 == 0) {
        cat("Processing position:", format(pos, big.mark=","), "\n")
      }
      
      # Process each sample at this position
      sample_results <- samples %>%
        map(function(sample_name) {
          result <- estimate_haplotypes_list_format(pos, sample_name, observed_euchromatic, 
                                                  founders, h_cutoff, "adaptive", 
                                                  NULL, chr, verbose = 0)
          return(result)
        }) %>%
        compact()  # Remove NULL results
      
      if (length(sample_results) == 0) {
        return(NULL)
      }
      
      # Combine results for all samples at this position
      all_samples <- unlist(map(sample_results, "sample"))
      all_groups <- map(sample_results, ~ .x$Groups[[1]][[1]])
      all_haps <- map(sample_results, ~ .x$Haps[[1]][[1]])
      all_err <- map(sample_results, ~ .x$Err[[1]][[1]])
      
      # Create the combined list format structure for this position
      tibble(
        CHROM = chr,
        pos = pos,
        sample = list(all_samples),
        Groups = list(all_groups),
        Haps = list(all_haps),
        Err = list(all_err),
        Names = list(rep(list(founders), length(all_samples)))
      )
    }) %>%
    filter(!is.null(pos))  # Remove NULL results
  
  # Save adaptive_h4 results
  output_file <- file.path(hap_list_results_dir, paste0("adaptive_h4_list_format_", chr, ".RDS"))
  saveRDS(list_format_data, output_file)
  cat("✓ Adaptive_h4 list format results saved:", nrow(list_format_data), "rows\n")
  cat("Output file:", output_file, "\n")
  
} else if (method == "smooth_h4") {
  # Run smooth_h4 with list format
  cat("Running smooth_h4 estimation...\n")
  
  # First, load the adaptive_h4 results to smooth over
  adaptive_file <- file.path(hap_list_results_dir, paste0("adaptive_h4_list_format_", chr, ".RDS"))
  if (!file.exists(adaptive_file)) {
    stop("Adaptive_h4 results not found. Run adaptive_h4 first: ", adaptive_file)
  }
  
  adaptive_data <- readRDS(adaptive_file)
  cat("✓ Loaded adaptive_h4 results:", nrow(adaptive_data), "positions\n")
  
  # Create smooth_h4 by averaging adaptive_h4 results
  # For smooth_h4, we use a sliding window approach to average nearby positions
  window_size <- 5  # Average over ±2 positions (5 total)
  
  smooth_data <- adaptive_data %>%
    mutate(
      # For smooth_h4, groups are always 1:8 (all founders distinguishable)
      Groups = list(rep(list(1:8), length(sample[[1]]))),
      
      # Average the haplotype frequencies using sliding window
      Haps = list(map(seq_along(sample[[1]]), function(i) {
        # Get haplotype estimates for this sample at nearby positions
        pos_idx <- which(adaptive_data$pos == pos)
        start_idx <- max(1, pos_idx - 2)
        end_idx <- min(nrow(adaptive_data), pos_idx + 2)
        
        # Get haplotype estimates for this sample across the window
        window_haps <- map(start_idx:end_idx, function(j) {
          adaptive_data$Haps[[j]][[1]][[i]]
        })
        
        # Average the haplotype frequencies
        valid_haps <- window_haps[map_lgl(window_haps, ~ !any(is.na(.x)))]
        
        if (length(valid_haps) > 0) {
          # Average valid estimates
          avg_haps <- reduce(valid_haps, `+`) / length(valid_haps)
          names(avg_haps) <- founders
          return(avg_haps)
        } else {
          # No valid estimates, return NAs
          return(set_names(rep(NA, length(founders)), founders))
        }
      })),
      
      # Average the error matrices using sliding window
      Err = list(map(seq_along(sample[[1]]), function(i) {
        # Get error matrices for this sample at nearby positions
        pos_idx <- which(adaptive_data$pos == pos)
        start_idx <- max(1, pos_idx - 2)
        end_idx <- min(nrow(adaptive_data), pos_idx + 2)
        
        # Get error matrices for this sample across the window
        window_errs <- map(start_idx:end_idx, function(j) {
          adaptive_data$Err[[j]][[1]][[i]]
        })
        
        # Average the error matrices
        valid_errs <- window_errs[map_lgl(window_errs, ~ !any(is.na(.x)))]
        
        if (length(valid_errs) > 0) {
          # Average valid error matrices
          avg_err <- reduce(valid_errs, `+`) / length(valid_errs)
          rownames(avg_err) <- founders
          colnames(avg_err) <- founders
          return(avg_err)
        } else {
          # No valid error matrices, return NAs
          return(matrix(NA, length(founders), length(founders), 
                       dimnames = list(founders, founders)))
        }
      }))
    )
  
  # Save smooth_h4 results
  output_file <- file.path(hap_list_results_dir, paste0("smooth_h4_list_format_", chr, ".RDS"))
  saveRDS(smooth_data, output_file)
  cat("✓ Smooth_h4 list format results saved:", nrow(smooth_data), "rows\n")
  cat("Output file:", output_file, "\n")
}

cat("\n=== HAPLOTYPE ESTIMATION COMPLETE ===\n")
cat("Results saved to:", hap_list_results_dir, "\n")
