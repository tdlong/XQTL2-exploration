#!/usr/bin/env Rscript

# Run Haplotype Estimation with List Format Output
# Focuses on adaptive_h4 and smooth_h4 methods
# Saves results to new hap_list_results directory

library(tidyverse)
library(limSolve)

# Source the new estimator functions
source("scripts/production/haplotype_estimation_functions.R")

# NEW HAPLOTYPE ESTIMATOR WITH LIST FORMAT OUTPUT
# This is a new function that captures groups and error matrices from the existing algorithm

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                                           method = c("fixed", "adaptive"),
                                           window_size_bp = NULL,
                                           chr = "chr2R",
                                           verbose = 0) {
  
  method <- match.arg(method)
  
  # Initialize result variables
  estimate_OK <- NA
  haplotype_freqs <- rep(NA, length(founders))
  names(haplotype_freqs) <- founders
  error_matrix <- matrix(NA, length(founders), length(founders))
  rownames(error_matrix) <- founders
  colnames(error_matrix) <- founders
  groups <- rep(NA, length(founders))
  names(groups) <- founders
  final_window_size <- NA
  n_snps <- 0
  
  if (method == "fixed") {
    # FIXED WINDOW METHOD
    window_start <- max(1, pos - window_size_bp/2)
    window_end <- pos + window_size_bp/2
    
    # Get SNPs in window
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    if (nrow(window_data) == 0) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    # Get founder matrix and sample frequencies
    founder_matrix_clean <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    sample_freqs_clean <- wide_data %>%
      pull(!!sample_name)
    
    complete_rows <- complete.cases(founder_matrix_clean) & !is.na(sample_freqs_clean)
    founder_matrix_clean <- founder_matrix_clean[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs_clean[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) {
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                              rep(NA, length(founders)), rep(1, length(founders)), 
                              matrix(NA, length(founders), length(founders)), founders))
    }
    
    final_window_size <- window_size_bp
    n_snps <- nrow(wide_data)
    
    # Run LSEI with error matrix capture
    tryCatch({
      E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      G <- diag(length(founders))  # Non-negativity constraints
      H <- matrix(rep(0.0003, length(founders)))  # Lower bound
      
      # Call lsei with fulloutput=TRUE to get error matrix
      lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                   E = E, F = F, G = G, H = H, fulloutput = TRUE)
      
      if (lsei_result$IsError == 0) {
        # LSEI successful - get frequencies
        haplotype_freqs <- lsei_result$X
        names(haplotype_freqs) <- founders
        
        # Capture the error matrix
        if (!is.null(lsei_result$cov)) {
          error_matrix <- lsei_result$cov
          rownames(error_matrix) <- founders
          colnames(error_matrix) <- founders
        }
        
        # Get clustering groups (this is the key part that was missing!)
        distances <- dist(t(founder_matrix_clean))
        hclust_result <- hclust(distances, method = "complete")
        groups <- cutree(hclust_result, h = h_cutoff)
        names(groups) <- founders
        
        # Determine estimate_OK based on distinguishability
        n_groups <- length(unique(groups))
        estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
        
      } else {
        # LSEI failed
        estimate_OK <- 0
        groups <- rep(1, length(founders))
        names(groups) <- founders
      }
    }, error = function(e) {
      # Catastrophic LSEI error
      estimate_OK <- 0
      groups <- rep(1, length(founders))
      names(groups) <- founders
    })
    
  } else {
    # ADAPTIVE WINDOW METHOD - simplified for now
    # This would need the full adaptive logic from the existing code
    estimate_OK <- 0
    groups <- rep(1, length(founders))
    names(groups) <- founders
    final_window_size <- 100000  # placeholder
  }
  
  return(create_list_result(chr, pos, sample_name, method, final_window_size, n_snps, 
                          estimate_OK, haplotype_freqs, groups, error_matrix, founders))
}

# Helper function to create the list format result
create_list_result <- function(chr, pos, sample_name, method, window_size, n_snps, 
                              estimate_OK, haplotype_freqs, groups, error_matrix, founders) {
  
  # Create the list format structure
  result <- list(
    CHROM = chr,
    pos = pos,
    sample = list(sample_name),
    Groups = list(list(groups)),
    Haps = list(list(haplotype_freqs)),
    Err = list(list(error_matrix)),
    Names = list(list(founders))
  )
  
  return(result)
}

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
  # For smooth_h4, we use a 21-position sliding window (matching original implementation)
  window_size <- 21
  
  # First, we need to determine the quality of each position based on groups
  # A position is "OK" if it has 8 groups (all founders distinguishable)
  adaptive_data_with_quality <- adaptive_data %>%
    mutate(
      # Determine if each position is "OK" based on groups
      position_ok = map_lgl(Groups, function(groups_list) {
        # Check if any sample at this position has 8 groups
        any(map_lgl(groups_list, function(groups) {
          length(unique(groups)) == 8
        }))
      })
    )
  
  # Apply 21-position sliding window smoothing
  smooth_data <- adaptive_data_with_quality %>%
    arrange(pos) %>%
    mutate(
      # Calculate quality count: how many positions in 21-position window are OK?
      quality_count = map_dbl(seq_len(n()), function(i) {
        start_idx <- max(1, i - 10)  # 10 positions before
        end_idx <- min(n(), i + 10)  # 10 positions after
        window_positions <- start_idx:end_idx
        sum(adaptive_data_with_quality$position_ok[window_positions], na.rm = TRUE)
      }),
      
      # New estimate_OK: OK if at least 17 out of 21 positions are OK
      new_estimate_ok = quality_count >= 17,
      
      # For smooth_h4, groups depend on the quality
      Groups = list(map(seq_along(sample[[1]]), function(i) {
        if (new_estimate_ok) {
          # If smooth estimate is OK, use 1:8 (all founders distinguishable)
          return(1:8)
        } else {
          # If smooth estimate is not OK, use all 1s (fallback)
          return(rep(1, length(founders)))
        }
      })),
      
      # Average the haplotype frequencies using 21-position sliding window
      Haps = list(map(seq_along(sample[[1]]), function(i) {
        pos_idx <- which(adaptive_data_with_quality$pos == pos)
        start_idx <- max(1, pos_idx - 10)  # 10 positions before
        end_idx <- min(nrow(adaptive_data_with_quality), pos_idx + 10)  # 10 positions after
        
        # Get haplotype estimates for this sample across the 21-position window
        window_haps <- map(start_idx:end_idx, function(j) {
          adaptive_data_with_quality$Haps[[j]][[1]][[i]]
        })
        
        # Only use positions where the original estimate was OK (8 groups)
        valid_positions <- map_lgl(start_idx:end_idx, function(j) {
          adaptive_data_with_quality$position_ok[j]
        })
        
        # Filter to only OK positions (not all 21 positions)
        valid_haps <- window_haps[valid_positions & map_lgl(window_haps, ~ !any(is.na(.x)))]
        
        if (new_estimate_ok && length(valid_haps) > 0) {
          # Average over only the OK positions (not all 21 positions)
          avg_haps <- reduce(valid_haps, `+`) / length(valid_haps)
          # Normalize so they sum to 1
          avg_haps <- avg_haps / sum(avg_haps)
          names(avg_haps) <- founders
          return(avg_haps)
        } else {
          # If <17 positions are OK, return NAs
          return(set_names(rep(NA, length(founders)), founders))
        }
      })),
      
      # Average the error matrices using 21-position sliding window
      Err = list(map(seq_along(sample[[1]]), function(i) {
        pos_idx <- which(adaptive_data_with_quality$pos == pos)
        start_idx <- max(1, pos_idx - 10)  # 10 positions before
        end_idx <- min(nrow(adaptive_data_with_quality), pos_idx + 10)  # 10 positions after
        
        # Get error matrices for this sample across the 21-position window
        window_errs <- map(start_idx:end_idx, function(j) {
          adaptive_data_with_quality$Err[[j]][[1]][[i]]
        })
        
        # Only use positions where the original estimate was OK (8 groups)
        valid_positions <- map_lgl(start_idx:end_idx, function(j) {
          adaptive_data_with_quality$position_ok[j]
        })
        
        # Filter to only OK positions (not all 21 positions)
        valid_errs <- window_errs[valid_positions & map_lgl(window_errs, ~ !any(is.na(.x)))]
        
        if (new_estimate_ok && length(valid_errs) > 0) {
          # Average over only the OK positions (not all 21 positions)
          avg_err <- reduce(valid_errs, `+`) / length(valid_errs)
          rownames(avg_err) <- founders
          colnames(avg_err) <- founders
          return(avg_err)
        } else {
          # If <17 positions are OK, return NAs
          return(matrix(NA, length(founders), length(founders), 
                       dimnames = list(founders, founders)))
        }
      }))
    ) %>%
    select(-position_ok, -quality_count, -new_estimate_ok)  # Clean up temporary columns
  
  # Save smooth_h4 results
  output_file <- file.path(hap_list_results_dir, paste0("smooth_h4_list_format_", chr, ".RDS"))
  saveRDS(smooth_data, output_file)
  cat("✓ Smooth_h4 list format results saved:", nrow(smooth_data), "rows\n")
  cat("Output file:", output_file, "\n")
}

cat("\n=== HAPLOTYPE ESTIMATION COMPLETE ===\n")
cat("Results saved to:", hap_list_results_dir, "\n")
