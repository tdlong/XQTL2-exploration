#!/usr/bin/env Rscript

# =============================================================================
# BENCHMARK: WIDE FORMAT vs PRODUCTION
# =============================================================================
# 
# This script runs both the production and wide format versions
# on the same data and compares their performance
# 
# Usage: Rscript benchmark_wide_vs_production.R <chr> <output_dir> <param_file>

suppressPackageStartupMessages({
  library(tidyverse)
  library(microbenchmark)
})

# Source the functions we need (avoid main execution blocks)
source("scripts/ErrMatrix/est_haps_wide.R")

# Load production functions manually to avoid main execution
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# Copy the production function directly to avoid sourcing issues
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  method <- match.arg(method)
  
  # ADAPTIVE WINDOW METHOD
  if (verbose >= 2) {
    cat(sprintf("=== ADAPTIVE WINDOW: h_cutoff=%g, pos=%s, sample=%s ===\n", 
                h_cutoff, format(pos, big.mark=","), sample_name))
  }
  
  window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
  if (verbose >= 2) {
    cat("Window sizes to test:", paste(window_sizes/1000, "kb"), "\n")
  }
  
  final_result <- NULL
  final_n_groups <- 0
  final_window_size <- window_sizes[1]
  final_wide_data <- NULL
  previous_n_groups <- 0
  
  # CONSTRAINT ACCUMULATION (core of adaptive algorithm)
  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL
  
  for (window_idx in seq_along(window_sizes)) {
    window_size <- window_sizes[window_idx]
    window_start <- max(1, pos - window_size/2)
    window_end <- pos + window_size/2
    
    if (verbose >= 2) {
      cat(sprintf("  Window %d: %s", window_idx, 
                  ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp"))))
    }
    
    # Get SNPs in window
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    if (nrow(window_data) == 0) next
    
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next
    
    # Get founder matrix and sample frequencies
    founder_matrix <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    sample_freqs <- wide_data %>%
      pull(!!sample_name)
    
    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) next
    
    # Always update final window info (this is the last window we actually tried)
    final_window_size <- window_size
    final_wide_data <- wide_data
    
    # Hierarchical clustering
    founder_dist <- dist(t(founder_matrix_clean))
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    
    if (verbose >= 2) {
      cat(sprintf(" - %d SNPs, %d groups", nrow(founder_matrix_clean), n_groups))
    }
    
    # Check if clustering improved
    if (window_idx > 1 && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat(" - ✗ No improvement\n")
      next  # No improvement, try larger window
    }
    
    previous_n_groups <- n_groups
    
    # Build constraint matrix with accumulated constraints
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add accumulated constraints from previous windows
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) {
        cat(sprintf("    ✓ Added %d accumulated constraints from previous windows\n", nrow(accumulated_constraints)))
      }
    } else {
      if (verbose >= 2) {
        cat("    • No accumulated constraints yet (first meaningful window)\n")
      }
    }
    
    # Run LSEI with constraints
    tryCatch({
      result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                              E = E, F = F, 
                              G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE)
      
      if (result$IsError == 0) {
        # LSEI successful - accumulate constraints for next window
        current_constraints <- NULL
        current_constraint_values <- NULL
        
        if (verbose >= 2) {
          cat(sprintf(" - ✓ LSEI success, %d groups", n_groups))
        }
        
        unique_clusters <- unique(groups)
        for (cluster_id in unique_clusters) {
          cluster_founders <- which(groups == cluster_id)
          if (length(cluster_founders) > 1) {
            # Create constraint row for this group
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            # Calculate the actual group frequency from lsei result
            group_freq <- sum(result$X[cluster_founders])
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, group_freq)
          } else {
            # Single founder: lock their exact frequency
            founder_freq <- result$X[cluster_founders]
            
            # Create constraint: this founder = their exact frequency
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, founder_freq)
          }
        }
        
        # Update accumulated constraints for next window
        if (!is.null(current_constraints)) {
          accumulated_constraints <- current_constraints
          accumulated_constraint_values <- current_constraint_values
          if (verbose >= 2) {
            cat(sprintf(" - %d constraints", nrow(current_constraints)))
          }
        } else {
          accumulated_constraints <- NULL
          accumulated_constraint_values <- NULL
        }
        
        # Store the result (will be overwritten as we expand)
        final_result <- result
        final_n_groups <- n_groups
        
        # Check if all founders are separated
        if (n_groups == length(founders)) {
          if (verbose >= 2) {
            cat(" - ✅ SUCCESS: All founders distinguished\n")
          }
          break  # Success! Stop expanding
        }
      }
    }, error = function(e) {
      # LSEI error, try larger window
      if (verbose >= 2) {
        cat(sprintf("LSEI error: %s\n", e$message))
      }
    })
  }
  
  # Apply the correct rules for output
  if (!is.null(final_result)) {
    # LSEI was successful - ALWAYS return the frequency estimates
    founder_frequencies <- final_result$X
    names(founder_frequencies) <- founders
    
    # Check if founders are distinguishable to set trust level
    if (final_n_groups == length(founders)) {
      estimate_OK <- 1  # Founders distinguishable - we can TRUST the frequencies
    } else {
      estimate_OK <- 0  # Founders NOT distinguishable - we CANNOT trust the frequencies
    }
  } else {
    # Either insufficient SNPs OR LSEI failed/didn't converge
    founder_frequencies <- rep(NA, length(founders))
    names(founder_frequencies) <- founders
    estimate_OK <- NA
  }
  
  # Return appropriate error matrix
  if (!is.null(final_result)) {
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      error_matrix <- final_result$covar
    } else {
      # LSEI didn't provide covariance - create NA matrix
      error_matrix <- matrix(NA, length(founders), length(founders))
      rownames(error_matrix) <- founders
      colnames(error_matrix) <- founders
    }
  } else {
    # Create NA matrix when no result
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders
    colnames(error_matrix) <- founders
  }
  
  return(list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders))
}

# Load wide format functions
source("scripts/ErrMatrix/wide_format_optimization.R")

if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript benchmark_wide_vs_production.R <chr> <output_dir> <param_file>")
  }
  
  chr <- args[1]
  output_dir <- args[2]
  param_file <- args[3]
  
  cat("=== BENCHMARKING WIDE FORMAT vs PRODUCTION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Load parameters
  source(param_file)
  
  # Load RefAlt data in wide format for comparison
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  cat("Loading data...\n")
  wide_data <- read_RefAlt_wide(refalt_file, founders)
  df3_wide <- wide_data$df3_wide
  founders <- wide_data$founders
  samples <- wide_data$samples
  
  # Convert to long format for production version
  df3_long <- df3_wide %>%
    dplyr::select(POS, dplyr::all_of(founders), dplyr::all_of(samples)) %>%
    tidyr::pivot_longer(cols = c(dplyr::all_of(founders), dplyr::all_of(samples)), 
                       names_to = "name", values_to = "freq") %>%
    dplyr::arrange(POS)
  
  # Define test positions (limited for benchmarking)
  euchromatin_boundaries <- list(
    chr2L = c(82455, 22011009),
    chr2R = c(5398184, 24684540),
    chr3L = c(158639, 22962476),
    chr3R = c(4552934, 31845060),
    chrX = c(277911, 22628490)
  )
  
  euchrom_start <- euchromatin_boundaries[[chr]][1]
  euchrom_end <- euchromatin_boundaries[[chr]][2]
  
  scan_start <- ceiling(euchrom_start / step) * step
  scan_end <- floor(euchrom_end / step) * step
  all_positions <- seq(scan_start, scan_end, by = step)
  
  # Limit to first 50 positions for benchmarking
  test_positions <- head(all_positions, 50)
  test_samples <- head(samples, 3)  # Test with 3 samples
  
  cat("Testing with", length(test_positions), "positions ×", length(test_samples), "samples\n")
  cat("Total operations:", length(test_positions) * length(test_samples), "\n\n")
  
  # Function to run production version
  run_production <- function() {
    results <- map_dfr(test_positions, function(test_pos) {
      map_dfr(test_samples, function(sample_name) {
        result <- estimate_haplotypes_list_format(
          pos = test_pos,
          sample_name = sample_name,
          df3 = df3_long,
          founders = founders,
          h_cutoff = 4,
          method = "adaptive",
          verbose = 0
        )
        
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = sample_name,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
    })
    return(results)
  }
  
  # Function to run wide format version
  run_wide_format <- function() {
    results <- map_dfr(test_positions, function(test_pos) {
      # Subset 500kb window for this position
      window_start <- max(1, test_pos - 250000)
      window_end <- min(max(df3_wide$POS), test_pos + 250000)
      
      window_data <- df3_wide %>%
        dplyr::filter(POS >= window_start & POS <= window_end)
      
      if (nrow(window_data) < 10) {
        # No data in window - return NA results for all samples
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = test_samples,
          Groups = map(test_samples, ~ rep(1, length(founders))),
          Haps = map(test_samples, ~ set_names(rep(NA, length(founders)), founders)),
          Err = map(test_samples, ~ matrix(NA, length(founders), length(founders))),
          Names = map(test_samples, ~ founders)
        ))
      }
      
      # Process all samples for this position
      sample_results <- map_dfr(test_samples, function(sample_name) {
        result <- est_haps_wide(
          pos = test_pos,
          sample_name = sample_name,
          df3_wide = window_data,
          founders = founders,
          h_cutoff = 4,
          verbose = 0
        )
        
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = sample_name,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
      
      return(sample_results)
    })
    return(results)
  }
  
  # Run benchmark
  cat("Running benchmark...\n")
  benchmark_results <- microbenchmark(
    production = run_production(),
    wide_format = run_wide_format(),
    times = 3,  # Run each version 3 times
    unit = "s"
  )
  
  # Display results
  cat("\n=== BENCHMARK RESULTS ===\n")
  print(benchmark_results)
  
  # Calculate speedup
  prod_times <- benchmark_results$time[benchmark_results$expr == "production"]
  wide_times <- benchmark_results$time[benchmark_results$expr == "wide_format"]
  
  prod_median <- median(prod_times) / 1e9  # Convert to seconds
  wide_median <- median(wide_times) / 1e9  # Convert to seconds
  
  speedup <- prod_median / wide_median
  
  cat("\n=== PERFORMANCE SUMMARY ===\n")
  cat("Production median time:", round(prod_median, 2), "seconds\n")
  cat("Wide format median time:", round(wide_median, 2), "seconds\n")
  cat("Speedup:", round(speedup, 2), "x\n")
  
  if (speedup > 1) {
    cat("✓ Wide format is", round(speedup, 2), "x FASTER\n")
  } else {
    cat("⚠ Wide format is", round(1/speedup, 2), "x SLOWER\n")
  }
  
  # Verify results are similar
  cat("\n=== VERIFYING RESULTS ===\n")
  prod_results <- run_production()
  wide_results <- run_wide_format()
  
  # Compare a few key metrics
  cat("Production results:", nrow(prod_results), "rows\n")
  cat("Wide format results:", nrow(wide_results), "rows\n")
  
  # Check if results have same structure
  if (nrow(prod_results) == nrow(wide_results)) {
    cat("✓ Same number of results\n")
  } else {
    cat("⚠ Different number of results\n")
  }
  
  cat("\n=== BENCHMARK COMPLETE ===\n")
}
