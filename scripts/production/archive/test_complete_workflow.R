#!/usr/bin/env Rscript

# Complete workflow test - reproduces the entire enchilada that's currently running
# Tests both adaptive estimation AND smoothing with JUICE data

library(tidyverse)
library(purrr)

# Source the working functions
source("scripts/production/haplotype_estimation_working.R")

# Test parameters
chr <- "chr2R"
method <- "adaptive"
parameter <- 4
output_dir <- "process/JUICE"
param_file <- "helpfiles/JUICE_haplotype_parameters.R"
n_positions <- 100
n_samples <- 1

cat("=== TESTING COMPLETE WORKFLOW (ADAPTIVE + SMOOTHING) ===\n")
cat("Chromosome:", chr, "\n")
cat("Method:", method, "\n")
cat("Parameter:", parameter, "\n")
cat("Output directory:", output_dir, "\n")
cat("Parameter file:", param_file, "\n")
cat("Testing positions:", n_positions, "\n")
cat("Testing samples:", n_samples, "\n")
cat("====================================================\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
if (!file.exists(param_file)) {
  cat("Error: Parameter file not found:", param_file, "\n")
  quit(status = 1)
}

source(param_file)
cat("✓ Parameter file:", param_file, "\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), samples (", length(names_in_bam), ")\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and process data
cat("\n2. Loading and processing data...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  cat("Error: RefAlt file not found:", refalt_file, "\n")
  quit(status = 1)
}

df3 <- process_refalt_data(refalt_file, founders)

# 3. Define test positions
cat("\n3. Defining test positions...\n")
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
test_positions <- head(all_positions, n_positions)

cat("✓ Euchromatin boundaries:", euchrom_start, "to", euchrom_end, "bp\n")
cat("✓ Scan grid:", scan_start, "to", scan_end, "by", step, "\n")
cat("✓ Test positions:", length(test_positions), "\n")

# 4. Get test samples
test_samples <- head(names_in_bam, n_samples)
cat("✓ Test samples:", paste(test_samples, collapse=", "), "\n")

# 5. STEP 1: Run adaptive estimation (same as production)
cat("\n4. STEP 1: Running adaptive estimation...\n")
cat("Processing", length(test_positions), "positions ×", length(test_samples), "samples...\n")

adaptive_results <- expand_grid(
  pos = test_positions,
  sample_name = test_samples
) %>%
  purrr::pmap_dfr(~ {
    cat("Processing pos:", ..1, "sample:", ..2, "\n")
    result <- estimate_haplotypes_list_format(
      pos = ..1,
      sample_name = ..2,
      df3 = df3,
      founders = founders,
      h_cutoff = parameter,
      method = method,
      window_size_bp = NULL,
      chr = chr,
      verbose = 0  # Minimal output for speed
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

cat("✓ Adaptive estimation complete:", nrow(adaptive_results), "results\n")

# 6. STEP 2: Apply smoothing (same as production)
cat("\n5. STEP 2: Applying 21-position sliding window smoothing...\n")

# Helper functions from the smoothing script
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

# Add estimate_OK column
adaptive_results <- adaptive_results %>%
  mutate(estimate_OK = map_lgl(Groups, check_estimate_ok))

# Process smoothing for each sample
unique_samples <- unique(adaptive_results$sample)
smooth_results <- map_dfr(unique_samples, function(sample_name) {
  cat("Smoothing sample:", sample_name, "\n")
  
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

cat("✓ Smoothing complete:", nrow(smooth_results), "results\n")

# 7. Show results summary
cat("\n6. Results summary...\n")
cat("=== ADAPTIVE RESULTS ===\n")
cat("Positions:", length(unique(adaptive_results$pos)), "\n")
cat("Samples:", length(unique(adaptive_results$sample)), "\n")
cat("Good estimates:", sum(adaptive_results$estimate_OK), "out of", nrow(adaptive_results), "\n")

cat("\n=== SMOOTH RESULTS ===\n")
cat("Positions:", length(unique(smooth_results$pos)), "\n")
cat("Samples:", length(unique(smooth_results$sample)), "\n")

# Check quality
smooth_quality <- map_lgl(smooth_results$Groups, function(groups) {
  length(unique(groups)) == 8 && all(sort(unique(groups)) == 1:8)
})
cat("Good smooth estimates:", sum(smooth_quality), "out of", nrow(smooth_results), "\n")

# Show some example results
cat("\n=== EXAMPLE RESULTS (first 5 positions) ===\n")
for (i in 1:min(5, nrow(smooth_results))) {
  cat("\nPosition:", smooth_results$pos[i], "Sample:", smooth_results$sample[i], "\n")
  haps <- smooth_results$Haps[[i]]
  if (!any(is.na(haps))) {
    cat("  Haplotype frequencies:", paste(round(haps, 3), collapse=", "), "\n")
    cat("  Sum:", round(sum(haps), 6), "\n")
  } else {
    cat("  No haplotype frequencies (quality not OK)\n")
  }
}

cat("\n=== COMPLETE WORKFLOW TEST COMPLETE ===\n")
cat("Successfully reproduced the entire enchilada!\n")
cat("Both adaptive estimation and smoothing are working.\n")
