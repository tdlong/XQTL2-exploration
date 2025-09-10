#!/usr/bin/env Rscript

# =============================================================================
# HAPLOTYPE ESTIMATION PIPELINE - FAST VECTORIZED VERSION
# =============================================================================
# 
# VECTORIZED pipeline for haplotype estimation - no loops!
# Processes entire chromosomes efficiently using dplyr group_by operations.
#
# Usage as script:
#   Rscript haplotype_estimation_pipeline_fast.R <chr> <method> <parameter> <output_dir> <param_file>
#   Rscript haplotype_estimation_pipeline_fast.R smooth <chr> <param_file> <output_dir>
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(readr)
  library(limSolve)
})

# =============================================================================
# VECTORIZED HAPLOTYPE ESTIMATION FUNCTIONS
# =============================================================================

#' Estimate haplotypes for ALL positions and samples at once (vectorized)
#' @param df3 Data frame with founder frequencies
#' @param founders Vector of founder names
#' @param h_cutoff Hierarchical clustering cutoff
#' @param debug Debug mode
#' @return Data frame with results for all position-sample combinations
estimate_haplotypes_vectorized <- function(df3, founders, h_cutoff, debug = FALSE) {
  
  cat("Running VECTORIZED haplotype estimation...\n")
  cat("Processing", nrow(df3), "position-sample combinations\n")
  
  # Group by position and process all samples for each position at once
  results <- df3 %>%
    group_by(POS) %>%
    group_modify(~ {
      pos <- .y$POS
      
      if (debug && pos %% 1000 == 0) {
        cat("Processing position", pos, "\n")
      }
      
      # Get all samples for this position
      pos_data <- .x
      
      # Create design matrix X (n_samples x n_founders) for this position
      X <- pos_data %>%
        select(all_of(founders)) %>%
        as.matrix()
      
      # Create response vector y (n_samples) for this position
      y <- pos_data$freq
      
      # Remove any rows with missing data
      complete_rows <- complete.cases(X) & !is.na(y)
      X <- X[complete_rows, , drop = FALSE]
      y <- y[complete_rows]
      sample_names <- pos_data$sample[complete_rows]
      
      if (nrow(X) == 0) {
        # Return empty results for this position
        return(tibble(
          sample = character(0),
          Groups = list(),
          Haps = list(),
          Err = list(),
          Names = list()
        ))
      }
      
      # Calculate distance matrix between founders (once per position)
      founder_dist <- dist(t(X), method = "euclidean")
      
      # Perform clustering (once per position)
      hc <- hclust(founder_dist, method = "complete")
      groups <- cutree(hc, h = h_cutoff)
      
      # Debug output for this position
      if (debug) {
        cat("Position", pos, "- Founder distance matrix:\n")
        dist_matrix <- as.matrix(founder_dist)
        rownames(dist_matrix) <- founders
        colnames(dist_matrix) <- founders
        print(round(dist_matrix^2, 0))
        cat("Groups:", paste(groups, collapse = " "), "\n\n")
      }
      
      # Solve the linear system for ALL samples at this position
      # Constraint 1: frequencies must sum to 1
      E <- matrix(1, nrow = 1, ncol = length(founders))
      f <- 1
      
      # Constraint 2: frequencies must be non-negative (>= 0.0003)
      G <- diag(length(founders))
      h <- rep(0.0003, length(founders))
      
      # Solve for each sample
      sample_results <- map(seq_len(nrow(X)), function(i) {
        result <- lsei(A = X[i, , drop = FALSE], B = y[i], E = E, F = f, G = G, H = h, fulloutput = TRUE)
        
        if (result$IsError) {
          return(list(
            Groups = groups,
            Haps = rep(NA, length(founders)),
            Err = matrix(NA, nrow = length(founders), ncol = length(founders)),
            Names = founders
          ))
        }
        
        # Extract the solution
        founder_frequencies <- result$X
        
        # Ensure frequencies sum to 1 and are non-negative
        founder_frequencies <- pmax(founder_frequencies, 0.0003)
        founder_frequencies <- founder_frequencies / sum(founder_frequencies)
        
        return(list(
          Groups = groups,
          Haps = founder_frequencies,
          Err = result$cov,
          Names = founders
        ))
      })
      
      # Convert to tibble
      tibble(
        sample = sample_names,
        Groups = map(sample_results, ~ .x$Groups),
        Haps = map(sample_results, ~ .x$Haps),
        Err = map(sample_results, ~ .x$Err),
        Names = map(sample_results, ~ .x$Names)
      )
    }) %>%
    ungroup() %>%
    mutate(pos = POS) %>%
    select(pos, sample, Groups, Haps, Err, Names)
  
  return(results)
}

#' Create smooth haplotype estimator from adaptive results
create_smooth_haplotype_estimator <- function(chr, param_file, output_dir) {
  cat("=== CREATING SMOOTH HAPLOTYPE ESTIMATOR ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Parameter file:", param_file, "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Source parameter file
  source(param_file)
  
  # Define file paths
  adaptive_file <- file.path(output_dir, "haplotype_results_list_format", "adapt_h4", paste0("R.haps.", chr, ".out.rds"))
  smooth_file <- file.path(output_dir, "haplotype_results_list_format", "smooth_h4", paste0("R.haps.", chr, ".out.rds"))
  
  # Create output directory
  smooth_dir <- dirname(smooth_file)
  if (!dir.exists(smooth_dir)) {
    dir.create(smooth_dir, recursive = TRUE)
  }
  
  # Load adaptive results
  if (!file.exists(adaptive_file)) {
    stop("Adaptive results file not found: ", adaptive_file)
  }
  
  cat("Loading adaptive results...\n")
  adaptive_data <- readRDS(adaptive_file)
  cat("✓ Loaded", nrow(adaptive_data), "positions\n")
  
  # Apply 21-position sliding window smoothing
  cat("Applying 21-position sliding window smoothing...\n")
  smooth_data <- adaptive_data %>%
    group_by(sample) %>%
    arrange(pos) %>%
    mutate(
      # Calculate quality for each position (how many of 21 positions have 8 groups)
      quality_count = map_dbl(seq_along(pos), ~ {
        start_idx <- max(1, .x - 10)
        end_idx <- min(n(), .x + 10)
        window_groups <- Groups[start_idx:end_idx]
        sum(map_lgl(window_groups, ~ length(unique(.x)) == 8))
      }),
      quality_ok = quality_count >= 17,
      
      # Calculate smooth estimates
      smooth_groups = map(seq_along(pos), ~ {
        if (quality_ok[.x]) {
          # Quality OK: use 8 distinguishable groups
          rep(1:8, length.out = length(founders))
        } else {
          # Quality NOT OK: cluster all founders together
          rep(1, length(founders))
        }
      }),
      
      smooth_haps = map(seq_along(pos), ~ {
        if (quality_ok[.x]) {
          # Quality OK: average over good positions in window
          start_idx <- max(1, .x - 10)
          end_idx <- min(n(), .x + 10)
          window_haps <- Haps[start_idx:end_idx]
          good_haps <- window_haps[quality_ok[start_idx:end_idx]]
          
          if (length(good_haps) > 0) {
            # Average the good positions
            avg_haps <- map_dbl(seq_along(founders), ~ {
              mean(map_dbl(good_haps, ~ .x[.y]), na.rm = TRUE)
            })
            # Normalize to sum to 1
            avg_haps / sum(avg_haps)
          } else {
            rep(NA, length(founders))
          }
        } else {
          # Quality NOT OK: return NAs
          rep(NA, length(founders))
        }
      }),
      
      smooth_err = map(seq_along(pos), ~ {
        if (quality_ok[.x]) {
          # Quality OK: average error matrices over good positions
          start_idx <- max(1, .x - 10)
          end_idx <- min(n(), .x + 10)
          window_err <- Err[start_idx:end_idx]
          good_err <- window_err[quality_ok[start_idx:end_idx]]
          
          if (length(good_err) > 0) {
            # Average the error matrices
            avg_err <- matrix(0, nrow = length(founders), ncol = length(founders))
            for (err_mat in good_err) {
              avg_err <- avg_err + err_mat
            }
            avg_err / length(good_err)
          } else {
            matrix(NA, nrow = length(founders), ncol = length(founders))
          }
        } else {
          # Quality NOT OK: return NAs
          matrix(NA, nrow = length(founders), ncol = length(founders))
        }
      })
    ) %>%
    ungroup() %>%
    select(pos, sample, Groups = smooth_groups, Haps = smooth_haps, Err = smooth_err, Names)
  
  # Save smooth results
  cat("Saving smooth results...\n")
  saveRDS(smooth_data, smooth_file)
  cat("✓ Saved smooth results to:", smooth_file, "\n")
  
  # Also save in the old format for compatibility
  old_format_file <- file.path(output_dir, "haplotype_results_list_format", paste0("smooth_h4_results_", chr, ".RDS"))
  saveRDS(smooth_data, old_format_file)
  cat("✓ Saved old format to:", old_format_file, "\n")
  
  cat("✓ Smooth haplotype estimator created successfully\n")
  return(smooth_data)
}

# =============================================================================
# MAIN PIPELINE FUNCTIONS
# =============================================================================

#' Run haplotype estimation for a single chromosome (FAST VERSION)
run_haplotype_estimation <- function(chr, method, parameter, output_dir, param_file, debug = FALSE) {
  cat("=== HAPLOTYPE ESTIMATION (FAST VECTORIZED) ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Method:", method, "\n")
  cat("Parameter:", parameter, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n\n")
  
  # Source parameter file
  source(param_file)
  
  # Define file paths
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  output_file <- file.path(output_dir, "haplotype_results_list_format", "adapt_h4", paste0("R.haps.", chr, ".out.rds"))
  
  # Create output directory
  output_dir_path <- dirname(output_file)
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  # Load RefAlt data
  cat("Loading RefAlt data...\n")
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  df <- read.table(refalt_file, header = TRUE)
  cat("✓ Loaded", nrow(df), "positions\n")
  
  # Process data
  cat("Processing data...\n")
  df3 <- df %>%
    mutate(CHROM = chr) %>%
    pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
    mutate(
      RefAlt = str_sub(lab, 1, 3),
      name = str_sub(lab, 5)
    ) %>%
    select(-lab) %>%
    pivot_wider(names_from = RefAlt, values_from = count) %>%
    mutate(
      freq = REF / (REF + ALT),
      N = REF + ALT
    ) %>%
    select(-REF, -ALT) %>%
    pivot_wider(names_from = name, values_from = c(freq, N), names_sep = "_") %>%
    pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "value") %>%
    separate(lab, into = c("type", "sample"), sep = "_") %>%
    pivot_wider(names_from = type, values_from = value) %>%
    filter(!is.na(freq)) %>%
    select(CHROM, POS, sample = sample, freq = freq, N = N)
  
  # Quality control: keep only positions where ALL founders are "fixed"
  founder_wide <- df3 %>%
    select(-N) %>%
    pivot_wider(names_from = sample, values_from = freq) %>%
    select(-CHROM, -POS)
  
  # Get founder names from parameter file
  founders <- names_in_bam
  
  # CRITICAL: Ensure column order matches founders order (pivot_wider can reorder lexically)
  founder_wide <- founder_wide[, c("POS", founders)]
  
  # Quality filter: keep only positions where all founders are "fixed" (<3% or >97%)
  quality_filter <- apply(founder_wide[, founders], 1, function(row) {
    all(row < 0.03 | row > 0.97)
  })
  
  founder_matrix_clean <- founder_wide[quality_filter, ]
  cat("✓ Quality-filtered positions:", nrow(founder_matrix_clean), "\n")
  
  # Get positions to scan
  scan_positions <- sort(unique(founder_matrix_clean$POS))
  cat("✓ Using", length(scan_positions), "positions for haplotype estimation\n\n")
  
  # Filter df3 to only include quality-filtered positions
  df3_clean <- df3 %>%
    filter(POS %in% scan_positions)
  
  cat("✓ Filtered data to", nrow(df3_clean), "position-sample combinations\n\n")
  
  # Run VECTORIZED haplotype estimation
  results <- estimate_haplotypes_vectorized(df3_clean, founders, h_cutoff, debug)
  
  # Save results
  cat("Saving results...\n")
  saveRDS(results, output_file)
  cat("✓ Results saved to:", output_file, "\n")
  cat("✓ Haplotype estimation completed successfully\n")
  
  return(results)
}

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

# Check if running as script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage:\n")
    cat("  Haplotype estimation: Rscript haplotype_estimation_pipeline_fast.R <chr> <method> <parameter> <output_dir> <param_file>\n")
    cat("  Smooth estimator:     Rscript haplotype_estimation_pipeline_fast.R smooth <chr> <param_file> <output_dir>\n")
    cat("  Help:                 Rscript haplotype_estimation_pipeline_fast.R help\n")
    quit(status = 0)
  }
  
  if (args[1] == "help") {
    cat("HAPLOTYPE ESTIMATION PIPELINE (FAST VERSION)\n")
    cat("===========================================\n\n")
    cat("This script provides VECTORIZED haplotype estimation - no loops!\n\n")
    cat("COMMANDS:\n")
    cat("1. Haplotype Estimation:\n")
    cat("   Rscript haplotype_estimation_pipeline_fast.R <chr> <method> <parameter> <output_dir> <param_file>\n")
    cat("   Example: Rscript haplotype_estimation_pipeline_fast.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R\n\n")
    cat("2. Smooth Estimator:\n")
    cat("   Rscript haplotype_estimation_pipeline_fast.R smooth <chr> <param_file> <output_dir>\n")
    cat("   Example: Rscript haplotype_estimation_pipeline_fast.R smooth chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2\n\n")
    cat("PERFORMANCE:\n")
    cat("  - Uses dplyr group_by operations instead of loops\n")
    cat("  - Processes all samples for each position at once\n")
    cat("  - Should be 10-100x faster than the original version\n")
    quit(status = 0)
  }
  
  if (args[1] == "smooth") {
    if (length(args) != 4) {
      stop("Usage: Rscript haplotype_estimation_pipeline_fast.R smooth <chr> <param_file> <output_dir>")
    }
    create_smooth_haplotype_estimator(args[2], args[3], args[4])
  } else {
    if (length(args) < 5 || length(args) > 6) {
      stop("Usage: Rscript haplotype_estimation_pipeline_fast.R <chr> <method> <parameter> <output_dir> <param_file> [debug]")
    }
    debug_mode <- if (length(args) >= 6) args[6] == "debug" else FALSE
    run_haplotype_estimation(args[1], args[2], as.numeric(args[3]), args[4], args[5], debug_mode)
  }
}
