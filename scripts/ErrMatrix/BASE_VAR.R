#!/usr/bin/env Rscript

# =============================================================================
# COMPLETE HAPLOTYPE WORKFLOW - ALL FUNCTIONS FROM 49H AGO
# =============================================================================
# 
# This file contains ALL functions needed to reproduce the complete workflow
# that the SLURM script runs. It combines:
# 1. Adaptive haplotype estimation (from run_haplotype_estimation_list_format.R)
# 2. Smoothing (from create_smooth_haplotype_estimator_list_format.R)
#
# ALL FUNCTIONS ARE EXACTLY COPIED FROM THE 49H AGO WORKING CODE
# NO MODIFICATIONS, NO FIXES, NO BASTARDIZATION
#
# =============================================================================
# FUNCTION DOCUMENTATION
# =============================================================================
#
# CORE HAPLOTYPE ESTIMATION:
# ---------------------------
# estimate_haplotypes_list_format(pos, sample_name, df3, founders, h_cutoff, ...)
#   PURPOSE: Wrapper function that calls est_haps_var (advanced variance/covariance estimation)
#   INPUT:  - pos: test position (bp) - passed as testing_position to est_haps_var
#           - sample_name: sample to estimate haplotypes for
#           - df3: long-format data (POS, name, freq)
#           - founders: vector of founder names
#           - h_cutoff: clustering threshold (typically 4)
#   OUTPUT: List with Groups, Haps, Err, Names (HARDWIRED format)
#   LOGIC:  - Calls est_haps_var with testing_position = pos
#           - est_haps_var uses genomic distance-based windowing (10kb-500kb)
#           - Advanced progressive V matrix construction for variance/covariance
#           - Pooled covariance estimation for grouped founders
#           - Constraint accumulation across window sizes
#
# DATA PROCESSING:
# ----------------
# process_refalt_data(refalt_file, founders)
#   PURPOSE: Load and process RefAlt.txt files into df3 format
#   INPUT:  - refalt_file: path to RefAlt.txt file
#           - founders: vector of founder names
#   OUTPUT: df3 tibble with columns: CHROM, POS, name, freq, N
#   LOGIC:  - Reads RefAlt.txt (wide format: CHROM, POS, F1_REF, F1_ALT, ...)
#           - Pivots to long format, calculates frequencies (REF/(REF+ALT))
#           - Applies quality filter: keeps only positions where ALL founders
#             are fixed (< 3% or > 97% frequency)
#
# SMOOTHING FUNCTIONS:
# --------------------
# check_estimate_ok(groups)
#   PURPOSE: Check if haplotype estimation was successful
#   INPUT:  - groups: vector of cluster assignments
#   OUTPUT: Boolean (TRUE if all 8 founders distinguishable)
#   LOGIC:  - Returns TRUE if exactly 8 unique groups (1:8)
#
# average_haps(haps_list, founders)
#   PURPOSE: Average haplotype frequencies across multiple positions
#   INPUT:  - haps_list: list of frequency vectors
#           - founders: vector of founder names
#   OUTPUT: Named vector of averaged frequencies
#   LOGIC:  - Averages all frequency vectors, normalizes to sum=1
#
# average_err(err_list, founders)
#   PURPOSE: Average error matrices across multiple positions
#   INPUT:  - err_list: list of error matrices
#           - founders: vector of founder names
#   OUTPUT: Averaged error matrix
#   LOGIC:  - Averages all error matrices element-wise
#
# MAIN WORKFLOW FUNCTIONS:
# ------------------------
# run_adaptive_estimation(chr, method, parameter, output_dir, param_file, ...)
#   PURPOSE: Run adaptive haplotype estimation for entire chromosome
#   INPUT:  - chr: chromosome name (chr2L, chr2R, etc.)
#           - method: "adaptive" (only method supported)
#           - parameter: h_cutoff value (typically 4)
#           - output_dir: where to save results
#           - param_file: R script with founders, names_in_bam, step
#   OUTPUT: Tibble with columns: CHROM, pos, sample, Groups, Haps, Err, Names
#   LOGIC:  - Loads parameters and RefAlt data
#           - Defines euchromatin boundaries and test positions (every 1kb)
#           - For each position×sample: calls estimate_haplotypes_list_format
#           - Saves results to adaptive_window_h4_results_<chr>.RDS
#
# run_smoothing(chr, param_file, output_dir, adaptive_results, ...)
#   PURPOSE: Apply 21-position sliding window smoothing to adaptive results
#   INPUT:  - chr: chromosome name
#           - param_file: R script with founders
#           - output_dir: where to save results
#           - adaptive_results: output from run_adaptive_estimation
#   OUTPUT: Tibble with smoothed results
#   LOGIC:  - For each sample: processes positions in order
#           - For each position: looks at 21-position window (±10 positions)
#           - Quality check: requires ≥17/21 positions with successful estimation
#           - If quality OK: averages haplotypes and errors from valid positions
#           - If quality poor: sets all founders to group 1, frequencies to NA
#           - Saves to smooth_h4_results_<chr>.RDS and reshaped formats
#
# WORKFLOW RELATIONSHIPS:
# -----------------------
# 1. process_refalt_data() → converts raw data to df3 format
# 2. run_adaptive_estimation() → calls estimate_haplotypes_list_format() for each position×sample
# 3. estimate_haplotypes_list_format() → calls est_haps_var() with advanced variance/covariance estimation
# 4. run_smoothing() → takes adaptive results, applies 21-position smoothing
# 5. Main execution → runs both steps in sequence
#
# DATA FLOW:
# ----------
# RefAlt.txt → process_refalt_data() → df3 (long format)
# df3 → estimate_haplotypes_list_format() → est_haps_var() → Groups, Haps, Err, Names
# Multiple results → run_adaptive_estimation() → adaptive_results tibble
# adaptive_results → run_smoothing() → smooth_results tibble
# Both results → saved as .RDS files for downstream analysis
#
# USAGE:
# ------
# Rscript scripts/production/complete_haplotype_workflow.R <chr> <method> <parameter> <output_dir> <param_file> [--nonverbose]
#
# EXAMPLES:
# Rscript scripts/production/complete_haplotype_workflow.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R
# Rscript scripts/production/complete_haplotype_workflow.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R --nonverbose
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# =============================================================================
# EXACT WORKING FUNCTION FROM 49H AGO - NO MODIFICATIONS
# =============================================================================

# Load est_haps_var from haplotype_error_workbench.R
source("scripts/ErrMatrix/haplotype_error_workbench.R")

# Alias for compatibility - est_haps_var uses testing_position instead of pos
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  # Call est_haps_var with testing_position parameter
  return(est_haps_var(
    testing_position = pos,
    sample_name = sample_name,
    df3 = df3,
    founders = founders,
    h_cutoff = h_cutoff,
    method = method,
    window_size_bp = window_size_bp,
    chr = chr,
    verbose = verbose
  ))
}

# =============================================================================
# DATA PROCESSING FUNCTIONS (EXACT FROM 49H AGO)
# =============================================================================

process_refalt_data <- function(refalt_file, founders) {
  # Load RefAlt data and process into df3 format - EXACT from working code
  cat("Loading RefAlt data from:", refalt_file, "\n")
  
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Process into df3 format (one row per sample per position) - EXACT from working code
  df3 <- refalt_data %>%
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
    select(-c("REF", "ALT")) %>%
    as_tibble()
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  founder_wide <- df3 %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  quality_filtered_positions <- founder_wide %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    ) %>%
    pull(POS)
  
  # Filter to quality positions and include sample data
  df3 <- df3 %>%
    filter(POS %in% quality_filtered_positions)
  
  cat("✓ Processed", nrow(df3), "rows for", length(unique(df3$name)), "samples\n")
  return(df3)
}

# =============================================================================
# SMOOTHING FUNCTIONS (EXACT FROM 49H AGO)
# =============================================================================

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

# =============================================================================
# MAIN WORKFLOW FUNCTIONS
# =============================================================================

run_adaptive_estimation <- function(chr, method, parameter, output_dir, param_file, debug = FALSE, verbose = TRUE, debug_level = 0) {
  # Step 1: Adaptive haplotype estimation - EXACT from working code
  #
  # CRITICAL: The output data frame structure is HARDWIRED and CANNOT be changed
  # Required columns: CHROM, pos, sample, Groups, Haps, Err, Names
  # - Groups: list of integer vectors (cluster assignments)
  # - Haps: list of named numeric vectors (founder frequencies)
  # - Err: list of matrices (error/covariance estimates)  
  # - Names: list of character vectors (founder names)
  # Any changes to this structure will break downstream analysis
  
  cat("=== RUNNING ADAPTIVE HAPLOTYPE ESTIMATION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Method:", method, "\n")
  cat("Parameter:", parameter, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n")
  cat("Debug mode:", debug, "\n")
  cat("Verbose mode:", verbose, "\n\n")
  
  # Load parameters
  source(param_file)
  
  # Load RefAlt data
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  df3 <- process_refalt_data(refalt_file, founders)
  
  # Define positions - EXACT from working code
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
  
  if (debug) {
    all_positions <- head(all_positions, 500)  # Limit to 500 positions for benchmarking
    # names_in_bam <- head(names_in_bam, 1)      # Run all 4 samples for benchmarking
  }
  
  total_operations <- length(all_positions) * length(names_in_bam)
  cat("Processing", length(all_positions), "positions ×", length(names_in_bam), "samples\n")
  cat("Total operations:", total_operations, "\n")
  
  if (!debug && total_operations > 100) {
    cat("This may take several hours to days. Progress will be shown every 100 operations.\n")
    cat("Estimated time per operation: 2-5 seconds (varies by data complexity)\n")
    cat("Estimated total time:", round(total_operations * 3.5 / 3600, 1), "hours\n\n")
  }
  
  # Run adaptive estimation with progress tracking
  if (debug || total_operations <= 100) {
    # Small dataset - no progress tracking needed
    adaptive_results <- expand_grid(
      pos = all_positions,
      sample_name = names_in_bam
    ) %>%
      purrr::pmap_dfr(~ {
        if (debug) cat("Processing pos:", ..1, "sample:", ..2, "\n")
        
        result <- estimate_haplotypes_list_format(
          pos = ..1,
          sample_name = ..2,
          df3 = df3,
          founders = founders,
          h_cutoff = parameter,
          method = method,
          window_size_bp = NULL,
          chr = chr,
          verbose = debug_level  # Use explicit debug level from command line
        )
        
        return(tibble(
          CHROM = chr,
          pos = ..1,
          sample = ..2,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
  } else {
    # Large dataset - show progress
    start_time <- Sys.time()
    operation_count <- 0
    
    adaptive_results <- expand_grid(
      pos = all_positions,
      sample_name = names_in_bam
    ) %>%
      purrr::pmap_dfr(~ {
        operation_count <<- operation_count + 1
        
        # Show progress every 100 operations
        if (operation_count %% 100 == 0 || operation_count == total_operations) {
          elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
          rate <- operation_count / elapsed
          remaining <- (total_operations - operation_count) / rate
          percent_done <- round(100 * operation_count / total_operations, 1)
          
          cat(sprintf("Progress: %d/%d (%.1f%%) | Rate: %.1f ops/sec | Elapsed: %.1f min | Remaining: %.1f min\n",
                     operation_count, total_operations, percent_done, rate, elapsed/60, remaining/60))
        }
        
        result <- estimate_haplotypes_list_format(
          pos = ..1,
          sample_name = ..2,
          df3 = df3,
          founders = founders,
          h_cutoff = parameter,
          method = method,
          window_size_bp = NULL,
          chr = chr,
          verbose = debug_level  # Use explicit debug level from command line
        )
        
        return(tibble(
          CHROM = chr,
          pos = ..1,
          sample = ..2,
          Groups = list(result$Groups),
          Haps = list(result$Haps),
          Err = list(result$Err),
          Names = list(result$Names)
        ))
      })
    
    total_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("\n✓ Adaptive estimation completed in %.1f minutes (%.1f ops/sec)\n", 
                total_time/60, total_operations/total_time))
  }
  
  # Save adaptive results
  list_results_dir <- file.path(output_dir, "haplotype_results_list_format")
  dir.create(list_results_dir, recursive = TRUE, showWarnings = FALSE)
  
  adaptive_file <- file.path(list_results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
  saveRDS(adaptive_results, adaptive_file)
  
  cat("✓ Adaptive estimation complete:", nrow(adaptive_results), "results\n")
  cat("✓ Saved to:", adaptive_file, "\n")
  
  return(adaptive_results)
}

run_smoothing <- function(chr, param_file, output_dir, adaptive_results, verbose = TRUE) {
  # Step 2: Apply 21-position sliding window smoothing - EXACT from working code
  #
  # CRITICAL: The output data frame structure is HARDWIRED and CANNOT be changed
  # Required columns: CHROM, pos, sample, Groups, Haps, Err, Names
  # - Groups: list of integer vectors (cluster assignments)
  # - Haps: list of named numeric vectors (founder frequencies)
  # - Err: list of matrices (error/covariance estimates)
  # - Names: list of character vectors (founder names)
  # Any changes to this structure will break downstream analysis
  
  cat("\n=== RUNNING SMOOTHING ===\n")
  
  # Load parameters to get founders
  source(param_file)
  
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
    
    # Only process positions that have full 21-position window
    valid_positions <- 11:(n_positions - 10)
    
    if (length(valid_positions) == 0) {
      return(tibble())
    }
    
    sample_smooth <- sample_data[valid_positions, ] %>%
      select(CHROM, pos, sample) %>%
      mutate(
        window_quality = map_dbl(valid_positions, function(i) {
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
  
  # Reshape to one row per position with samples as lists
  smooth_data_reshaped <- smooth_results %>%
    dplyr::arrange(CHROM, pos, sample) %>%
    dplyr::group_by(CHROM, pos) %>%
    dplyr::summarise(
      sample = list(sample),
      Groups = list(Groups),
      Haps   = list(Haps),
      Err    = list(Err),
      Names  = list(Names),
      .groups = "drop"
    )
  
  # Also reshape the original adaptive_h4 data to the same format
  adaptive_data_reshaped <- adaptive_results %>%
    dplyr::arrange(CHROM, pos, sample) %>%
    dplyr::group_by(CHROM, pos) %>%
    dplyr::summarise(
      sample = list(sample),
      Groups = list(Groups),
      Haps   = list(Haps),
      Err    = list(Err),
      Names  = list(Names),
      .groups = "drop"
    )
  
  # Save results
  # CRITICAL: Output file formats and locations are HARDWIRED and CANNOT be changed
  # The downstream pipeline expects these exact file names and structures:
  # - adaptive_window_h4_results_<chr>.RDS (original format)
  # - smooth_h4_results_<chr>.RDS (original format)  
  # - smooth_h4/R.haps.<chr>.out.rds (reshaped format)
  # - adapt_h4/R.haps.<chr>.out.rds (reshaped format)
  # Any changes to file names or structures will break the pipeline
  list_results_dir <- file.path(output_dir, "haplotype_results_list_format")
  smooth_dir <- file.path(list_results_dir, "smooth_h4")
  adapt_dir  <- file.path(list_results_dir, "adapt_h4")
  dir.create(smooth_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(adapt_dir,  recursive = TRUE, showWarnings = FALSE)
  
  # Save smooth_h4 in both formats
  smooth_original_file <- file.path(list_results_dir, paste0("smooth_h4_results_", chr, ".RDS"))
  smooth_reshaped_file <- file.path(smooth_dir, paste0("R.haps.", chr, ".out.rds"))
  
  # Save adaptive_h4 in reshaped format
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

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: Rscript complete_haplotype_workflow.R <chr> <method> <parameter> <output_dir> <param_file> [--nonverbose] [--debug-level-1] [--debug-level-2] [--debug-level-3]")
  }
  
  chr <- args[1]
  method <- args[2]
  parameter <- as.numeric(args[3])
  output_dir <- args[4]
  param_file <- args[5]
  debug <- "--debug" %in% args
  verbose <- !("--nonverbose" %in% args)
  
  # Parse debug level from command line
  debug_level <- 0
  if ("--debug-level-1" %in% args) debug_level <- 1
  if ("--debug-level-2" %in% args) debug_level <- 2
  if ("--debug-level-3" %in% args) debug_level <- 3
  
  # Run the complete workflow
  if (debug) {
    cat("=== COMPLETE HAPLOTYPE WORKFLOW (DEBUG MODE) ===\n")
    cat("Limited to 500 positions × 4 samples for benchmarking\n")
  } else {
    cat("=== COMPLETE HAPLOTYPE WORKFLOW ===\n")
    cat("Processing all positions and samples\n")
  }
  
  if (verbose) {
    cat("Verbose output enabled\n")
  } else {
    cat("Minimal output mode\n")
  }
  cat("Debug level:", debug_level, "\n\n")
  
  # Step 1: Adaptive estimation ONLY (no smoothing) - WITH TIMING
  cat("Starting adaptive estimation timing...\n")
  start_time <- Sys.time()
  
  adaptive_results <- run_adaptive_estimation(chr, method, parameter, output_dir, param_file, debug, verbose, debug_level)
  
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  cat("\n=== TIMING RESULTS ===\n")
  cat("Total time:", round(as.numeric(total_time, units = "secs"), 2), "seconds\n")
  cat("Total time:", round(as.numeric(total_time, units = "mins"), 2), "minutes\n")
  if (nrow(adaptive_results) > 0) {
    cat("Function calls:", nrow(adaptive_results), "\n")
    cat("Time per call:", round(as.numeric(total_time, units = "secs") / nrow(adaptive_results), 4), "seconds\n")
  }
  cat("========================\n\n")
  
  # Step 2: Smoothing - SKIPPED
  # smooth_results <- run_smoothing(chr, param_file, output_dir, adaptive_results, verbose)
  
  cat("\n=== WORKFLOW COMPLETE ===\n")
  cat("✓ Adaptive estimation completed successfully\n")
  cat("✓ Smoothing skipped (as requested)\n")
  cat("✓ Output files created in production format\n")
}
