#!/usr/bin/env Rscript

# Run only the smoothing step on existing adaptive results
# This is much faster than rerunning the entire pipeline

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_smoothing_only.R <chr> <output_dir> <param_file>")
}

chr <- args[1]
output_dir <- args[2]
param_file <- args[3]

cat("=== RUNNING SMOOTHING ONLY ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n\n")

# Load the existing adaptive results
list_results_dir <- file.path(output_dir, "haplotype_results_list_format")
adaptive_file <- file.path(list_results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))

if (!file.exists(adaptive_file)) {
  stop("Adaptive results file not found: ", adaptive_file)
}

cat("Loading adaptive results from:", basename(adaptive_file), "\n")
adaptive_results <- readRDS(adaptive_file)
cat("Loaded", nrow(adaptive_results), "adaptive results\n\n")

# Load parameters
source(param_file)

# Define the smoothing functions (copied from BASE_VAR_WIDE.R)
check_estimate_ok <- function(groups) {
  if (is.null(groups) || any(is.na(groups))) return(FALSE)
  length(unique(groups)) == 8 && all(sort(unique(groups)) == 1:8)
}

average_haps <- function(haps_list, founders) {
  if (length(haps_list) == 0) {
    return(set_names(rep(NA, length(founders)), founders))
  }
  avg_haps <- reduce(haps_list, `+`) / length(haps_list)
  
  if (sum(avg_haps) > 1e-10) {
    avg_haps <- avg_haps / sum(avg_haps)
  } else {
    avg_haps <- rep(1/length(founders), length(founders))
  }
  
  names(avg_haps) <- founders
  return(avg_haps)
}

average_err <- function(err_list, founders) {
  if (length(err_list) == 0) {
    na_matrix <- matrix(NA, length(founders), length(founders))
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  first_err <- err_list[[1]]
  if (is.null(first_err) || !is.matrix(first_err)) {
    na_matrix <- matrix(NA, length(founders), length(founders))
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  n_founders <- length(founders)
  if (nrow(first_err) != n_founders || ncol(first_err) != n_founders) {
    na_matrix <- matrix(NA, n_founders, n_founders)
    rownames(na_matrix) <- founders
    colnames(na_matrix) <- founders
    return(na_matrix)
  }
  
  avg_err <- reduce(err_list, `+`) / length(err_list)
  rownames(avg_err) <- founders
  colnames(avg_err) <- founders
  return(avg_err)
}

# The corrected smoothing function
run_smoothing <- function(chr, param_file, output_dir, adaptive_results, verbose = TRUE) {
  cat("\n=== RUNNING SMOOTHING ===\n")

  # Add estimate_OK column
  adaptive_results <- adaptive_results %>%
    mutate(estimate_OK = map_lgl(Groups, check_estimate_ok))

  # Process smoothing for each sample
  unique_samples <- unique(adaptive_results$sample)
  total_samples <- length(unique_samples)

  if (total_samples > 1 && verbose) {
    cat("Smoothing", total_samples, "samples...\n")
  }

  smooth_results <- map_dfr(seq_along(unique_samples), function(sample_idx) {
    sample_name <- unique_samples[sample_idx]

    if (total_samples > 1 && verbose) {
      cat(sprintf("Smoothing sample %d/%d: %s\n", sample_idx, total_samples, sample_name))
    }

    sample_data <- adaptive_results %>% 
      filter(sample == sample_name) %>%
      arrange(pos)

    n_positions <- nrow(sample_data)

    # Only process positions that have full 21-position window (±10 from center)
    # Position i can be smoothed if we have positions i-10 through i+10
    # So valid positions are 11 through (n_positions - 10)
    valid_positions <- 11:(n_positions - 10)
    
    if (length(valid_positions) == 0) {
      return(tibble())
    }
    
    sample_smooth <- sample_data[valid_positions, ] %>%
      select(CHROM, pos, sample) %>%
      mutate(
        window_quality = map_dbl(valid_positions, function(i) {
          # For 21-position window (±10 from center), we need positions i-10 to i+10
          start_idx <- i - 10
          end_idx <- i + 10
          window_indices <- start_idx:end_idx
          sum(sample_data$estimate_OK[window_indices], na.rm = TRUE)
        }),
        
        quality_ok = window_quality >= 17,
        
        Groups = map2(valid_positions, quality_ok, function(i, quality) {
          if (quality) {
            return(1:8)
          } else {
            return(rep(1, length(founders)))
          }
        }),
        
        Haps = map2(valid_positions, quality_ok, function(i, quality) {
          if (!quality) {
            return(set_names(rep(NA, length(founders)), founders))
          }
          
          # Use full 21-position window
          start_idx <- i - 10
          end_idx <- i + 10
          window_indices <- start_idx:end_idx
          
          valid_haps <- sample_data$Haps[window_indices][sample_data$estimate_OK[window_indices]]
          valid_haps <- valid_haps[map_lgl(valid_haps, ~ !any(is.na(.x)))]
          
          return(average_haps(valid_haps, founders))
        }),
        
        Err = map2(valid_positions, quality_ok, function(i, quality) {
          if (!quality) {
            na_matrix <- matrix(NA, length(founders), length(founders))
            rownames(na_matrix) <- founders
            colnames(na_matrix) <- founders
            return(na_matrix)
          }
          
          # Use full 21-position window
          start_idx <- i - 10
          end_idx <- i + 10
          window_indices <- start_idx:end_idx
          
          valid_errs <- sample_data$Err[window_indices][sample_data$estimate_OK[window_indices]]
          valid_errs <- valid_errs[map_lgl(valid_errs, ~ !any(is.na(.x)))]
          
          return(average_err(valid_errs, founders))
        }),
        
        Names = map(valid_positions, ~ founders)
      ) %>%
      select(CHROM, pos, sample, Groups, Haps, Err, Names)

    return(sample_smooth)
  })

  # Reshape to one row per position with samples as lists (preserve alignment)
  smooth_data_reshaped <- smooth_results %>%
    dplyr::arrange(CHROM, pos, sample) %>%
    dplyr::group_by(CHROM, pos) %>%
    dplyr::summarise({
      ord <- order(sample)
      tibble(
        sample = list(sample[ord]),
        Groups = list(Groups[ord]),
        Haps   = list(Haps[ord]),
        Err    = list(Err[ord]),
        Names  = list(Names[ord])
      )
    }, .groups = "drop")

  # Also reshape the original adaptive_h4 data to the same format (preserve alignment)
  adaptive_data_reshaped <- adaptive_results %>%
    dplyr::arrange(CHROM, pos, sample) %>%
    dplyr::group_by(CHROM, pos) %>%
    dplyr::summarise({
      ord <- order(sample)
      tibble(
        sample = list(sample[ord]),
        Groups = list(Groups[ord]),
        Haps   = list(Haps[ord]),
        Err    = list(Err[ord]),
        Names  = list(Names[ord])
      )
    }, .groups = "drop")

  # Save results
  list_results_dir <- file.path(output_dir, "haplotype_results_list_format")
  smooth_dir <- file.path(list_results_dir, "smooth_h4")
  adapt_dir  <- file.path(list_results_dir, "adapt_h4")
  dir.create(smooth_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(adapt_dir,  recursive = TRUE, showWarnings = FALSE)

  smooth_original_file <- file.path(list_results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
  smooth_reshaped_file <- file.path(smooth_dir, paste0("R.haps.", chr, ".out.rds"))
  adaptive_reshaped_file <- file.path(adapt_dir, paste0("R.haps.", chr, ".out.rds"))

  saveRDS(smooth_results, smooth_original_file)
  saveRDS(smooth_data_reshaped, smooth_reshaped_file)
  saveRDS(adaptive_data_reshaped, adaptive_reshaped_file)

  cat("✓ Smoothing complete:", nrow(smooth_results), "results\n")
  cat("✓ Saved smooth_h4 (original) to:", basename(smooth_original_file), "\n")
  cat("✓ Saved smooth_h4 (reshaped) to:", basename(smooth_reshaped_file), "\n")
  cat("✓ Saved adaptive_h4 (reshaped) to:", basename(adaptive_reshaped_file), "\n")

  return(smooth_results)
}

# Run smoothing
cat("Running smoothing...\n")
smooth_results <- run_smoothing(chr, param_file, output_dir, adaptive_results, verbose = TRUE)

cat("\n=== SMOOTHING COMPLETE ===\n")
cat("Smoothed results:", nrow(smooth_results), "rows\n")

# Check the output files
cat("\n=== CHECKING OUTPUT FILES ===\n")

# Check smooth_h4_results file
smooth_original_file <- file.path(list_results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
if (file.exists(smooth_original_file)) {
  smooth_original <- readRDS(smooth_original_file)
  cat("smooth_h4_results_", chr, ".RDS: ", nrow(smooth_original), " rows\n", sep="")
} else {
  cat("smooth_h4_results_", chr, ".RDS: NOT FOUND\n", sep="")
}

# Check smooth_h4 reshaped file
smooth_reshaped_file <- file.path(list_results_dir, "smooth_h4", paste0("R.haps.", chr, ".out.rds"))
if (file.exists(smooth_reshaped_file)) {
  smooth_reshaped <- readRDS(smooth_reshaped_file)
  cat("smooth_h4/R.haps.", chr, ".out.rds: ", nrow(smooth_reshaped), " rows\n", sep="")
} else {
  cat("smooth_h4/R.haps.", chr, ".out.rds: NOT FOUND\n", sep="")
}

# Check adapt_h4 reshaped file
adaptive_reshaped_file <- file.path(list_results_dir, "adapt_h4", paste0("R.haps.", chr, ".out.rds"))
if (file.exists(adaptive_reshaped_file)) {
  adaptive_reshaped <- readRDS(adaptive_reshaped_file)
  cat("adapt_h4/R.haps.", chr, ".out.rds: ", nrow(adaptive_reshaped), " rows\n", sep="")
} else {
  cat("adapt_h4/R.haps.", chr, ".out.rds: NOT FOUND\n", sep="")
}

# Compare with original
cat("\n=== POSITION COMPARISON ===\n")
orig_positions <- length(unique(adaptive_results$pos))
smooth_positions <- length(unique(smooth_results$pos))
expected_smooth <- orig_positions - 20

cat("Original positions:", orig_positions, "\n")
cat("Smoothed positions:", smooth_positions, "\n")
cat("Expected smoothed positions:", expected_smooth, "\n")
cat("Difference:", orig_positions - smooth_positions, "(should be 20)\n")

if (smooth_positions == expected_smooth) {
  cat("✅ CORRECT: Smoothed data has exactly 20 fewer positions\n")
} else {
  cat("❌ ERROR: Smoothed data has wrong number of positions\n")
}

cat("\n=== DONE ===\n")