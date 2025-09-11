#!/usr/bin/env Rscript

# =============================================================================
# WIDE FORMAT OPTIMIZATION - PRODUCTION COPY FOR OPTIMIZATION
# =============================================================================
# 
# This file contains copies of production functions needed for wide-format optimization:
# 1. process_refalt_data() - modified to work with wide format
# 2. run_adaptive_estimation() - modified to use wide format and batch processing
# 3. estimate_haplotypes_list_format() - placeholder for new optimized version
#
# GOAL: Optimize the production pipeline by:
# - Reading entire chromosome in wide format
# - Converting counts → frequencies + QC filtering (wide format)
# - For each 10kb test position: subset 500kb window, map over all samples
# - Eliminate repeated pivot_wider and arrange operations
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# =============================================================================
# DATA PROCESSING FUNCTIONS (MODIFIED FOR WIDE FORMAT)
# =============================================================================

read_RefAlt_wide <- function(refalt_file, founders) {
  # Load RefAlt data and process into WIDE format for optimization
  # INPUT:  - refalt_file: path to RefAlt.txt file
  #         - founders: vector of founder names
  # OUTPUT: List with:
  #         - df3_wide: wide format data (POS, F1, F2, ..., sample1, sample2, ...)
  #         - founders: founder names
  #         - samples: sample names
  #         - quality_positions: positions that passed QC
  
  cat("Loading RefAlt data from:", refalt_file, "\n")
  
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Process into wide format directly (no pivot_longer needed)
  # Convert REF/ALT counts to frequencies
  founder_cols <- paste0(founders, "_REF")
  alt_cols <- paste0(founders, "_ALT")
  
  # Calculate frequencies for founders
  founder_freqs <- refalt_data[founder_cols] / (refalt_data[founder_cols] + refalt_data[alt_cols])
  names(founder_freqs) <- founders
  
  # Get sample columns (everything that's not CHROM, POS, or founder REF/ALT)
  all_cols <- names(refalt_data)
  sample_cols <- all_cols[!all_cols %in% c("CHROM", "POS", founder_cols, alt_cols)]
  
  # Calculate frequencies for samples
  sample_freqs <- refalt_data[sample_cols] / (refalt_data[sample_cols] + refalt_data[gsub("_REF", "_ALT", sample_cols)])
  
  # Create wide format dataframe
  df3_wide <- data.frame(
    CHROM = refalt_data$CHROM,
    POS = refalt_data$POS,
    founder_freqs,
    sample_freqs
  )
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  quality_mask <- apply(founder_freqs, 1, function(row) {
    all(is.na(row) | row < 0.03 | row > 0.97)
  })
  
  # Filter to quality positions
  df3_wide <- df3_wide[quality_mask, ]
  
  # Extract sample names (remove _REF suffix)
  sample_names <- gsub("_REF", "", sample_cols)
  
  cat("✓ Processed", nrow(df3_wide), "positions for", length(founders), "founders and", length(sample_names), "samples\n")
  cat("✓ Quality filter kept", nrow(df3_wide), "positions\n")
  
  return(list(
    df3_wide = df3_wide,
    founders = founders,
    samples = sample_names
  ))
}

# =============================================================================
# HAPLOTYPE ESTIMATION (USE EXISTING FROM WORKBENCH)
# =============================================================================

# We'll use the existing haplotype estimator from the workbench
# This will be aliased as est_haps_wide in the main execution section

# =============================================================================
# OPTIMIZED ADAPTIVE ESTIMATION (MODIFIED FOR WIDE FORMAT)
# =============================================================================

run_adapt_h4_wide <- function(chr, method, parameter, output_dir, param_file, debug = FALSE, verbose = TRUE, debug_level = 0) {
  # Modified version of run_adaptive_estimation for wide format optimization
  # 
  # OPTIMIZATION STRATEGY:
  # 1. Load entire chromosome in wide format (no pivot_wider)
  # 2. For each 10kb test position:
  #    - Subset 500kb window from wide matrix
  #    - Map haplotype estimation over all samples for that position
  #    - Store results for all samples at that position
  # 3. Eliminate repeated data reshaping operations
  
  cat("=== RUNNING OPTIMIZED ADAPTIVE HAPLOTYPE ESTIMATION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Method:", method, "\n")
  cat("Parameter:", parameter, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n")
  cat("Debug mode:", debug, "\n")
  cat("Verbose mode:", verbose, "\n\n")
  
  # Load parameters
  source(param_file)
  
  # Load RefAlt data in wide format
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  wide_data <- read_RefAlt_wide(refalt_file, founders)
  df3_wide <- wide_data$df3_wide
  founders <- wide_data$founders
  samples <- wide_data$samples
  
  # Define positions (same as production)
  euchromatin_boundaries <- list(
    chr2L = c(82455, 22011009),
    chr2R = c(5398184, 24684540),
    chr3L = c(158639, 22962476),
    chr3R = c(4552934, 31845060),
    chrX = c(277911, 22628490)
  )
  
  euchrom_start <- euchromatin_boundaries[[chr]][1]
  euchrom_end <- euchromatin_boundaries[[chr]][2]
  
  # Get actual max position from data (in case it's less than euchrom_end)
  max_pos <- max(df3_wide$POS)
  
  scan_start <- ceiling(euchrom_start / step) * step
  scan_end <- floor(euchrom_end / step) * step
  all_positions <- seq(scan_start, scan_end, by = step)
  
  if (debug) {
    all_positions <- head(all_positions, 100)  # Limit to 100 positions for debugging
    samples <- head(samples, 1)                # Limit to 1 sample for debugging
  }
  
  total_operations <- length(all_positions) * length(samples)
  cat("Processing", length(all_positions), "positions ×", length(samples), "samples\n")
  cat("Total operations:", total_operations, "\n")
  
  if (!debug && total_operations > 100) {
    cat("This may take several hours to days. Progress will be shown every 100 operations.\n")
    cat("Estimated time per operation: 2-5 seconds (varies by data complexity)\n")
    cat("Estimated total time:", round(total_operations * 3.5 / 3600, 1), "hours\n\n")
  }
  
  # OPTIMIZED PROCESSING:
  # For each test position, subset 500kb window, then map over all samples
  if (debug || total_operations <= 100) {
    # Small dataset - no progress tracking needed
    adaptive_results <- map_dfr(all_positions, function(test_pos) {
      if (debug) cat("Processing position:", test_pos, "\n")
      
      # Subset 500kb window for this position (centered)
      window_start <- max(1, test_pos - 250000)
      window_end <- min(max_pos, test_pos + 250000)
      
      window_data <- df3_wide[df3_wide$POS >= window_start & df3_wide$POS <= window_end, ]
      
      if (nrow(window_data) < 10) {
        # No data in window - return NA results for all samples
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = samples,
          Groups = map(samples, ~ rep(1, length(founders))),
          Haps = map(samples, ~ set_names(rep(NA, length(founders)), founders)),
          Err = map(samples, ~ matrix(NA, length(founders), length(founders))),
          Names = map(samples, ~ founders)
        ))
      }
      
      # Process all samples for this position using the current haplotype estimator
      sample_results <- map_dfr(samples, function(sample_name) {
        result <- est_haps_wide(
          pos = test_pos,
          sample_name = sample_name,
          df3 = window_data,  # Note: workbench expects df3 format
          founders = founders,
          h_cutoff = parameter
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
  } else {
    # Large dataset - show progress
    start_time <- Sys.time()
    operation_count <- 0
    
    adaptive_results <- map_dfr(all_positions, function(test_pos) {
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
      
      # Subset 500kb window for this position (centered)
      window_start <- max(1, test_pos - 250000)
      window_end <- min(max_pos, test_pos + 250000)
      
      window_data <- df3_wide[df3_wide$POS >= window_start & df3_wide$POS <= window_end, ]
      
      if (nrow(window_data) < 10) {
        # No data in window - return NA results for all samples
        return(tibble(
          CHROM = chr,
          pos = test_pos,
          sample = samples,
          Groups = map(samples, ~ rep(1, length(founders))),
          Haps = map(samples, ~ set_names(rep(NA, length(founders)), founders)),
          Err = map(samples, ~ matrix(NA, length(founders), length(founders))),
          Names = map(samples, ~ founders)
        ))
      }
      
      # Process all samples for this position using the current haplotype estimator
      sample_results <- map_dfr(samples, function(sample_name) {
        result <- est_haps_wide(
          pos = test_pos,
          sample_name = sample_name,
          df3 = window_data,  # Note: workbench expects df3 format
          founders = founders,
          h_cutoff = parameter
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
    
    total_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("\n✓ Optimized estimation completed in %.1f minutes (%.1f ops/sec)\n", 
                total_time/60, total_operations/total_time))
  }
  
  # Save results in a SAFE directory (won't overwrite production data)
  list_results_dir <- file.path(output_dir, "haplotype_results_list_format_wide_optimized")
  dir.create(list_results_dir, recursive = TRUE, showWarnings = FALSE)
  
  adaptive_file <- file.path(list_results_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
  saveRDS(adaptive_results, adaptive_file)
  
  cat("✓ Optimized estimation complete:", nrow(adaptive_results), "results\n")
  cat("✓ Saved to SAFE directory:", adaptive_file, "\n")
  cat("✓ This will NOT overwrite existing production data\n")
  
  return(adaptive_results)
}

# =============================================================================
# MAIN EXECUTION (FOR TESTING) - MOVED TO SEPARATE FILE
# =============================================================================

# Main execution block moved to test_wide_format.R to avoid sourcing issues
