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
# USAGE:
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

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  # CRITICAL: This function's return format is HARDWIRED and CANNOT be changed
  # The downstream pipeline expects exactly this structure:
  # - Groups: integer vector of cluster assignments
  # - Haps: named numeric vector of founder frequencies  
  # - Err: matrix of error/covariance estimates
  # - Names: character vector of founder names
  # Any changes to this format will break the entire pipeline
  
  # Debug output only in verbose mode
  if (verbose >= 1) {
    cat("DEBUG: Using estimate_haplotypes_list_format function\n")
  }
  
  method <- match.arg(method)
  
  if (verbose >= 1) {
    cat(sprintf("Processing pos: %s, sample: %s, method: %s\n", 
                format(pos, big.mark=","), sample_name, method))
  }
    # ADAPTIVE WINDOW METHOD
    if (verbose >= 2) {
      cat(sprintf("=== ADAPTIVE WINDOW METHOD: h_cutoff = %g ===\n", h_cutoff))
      cat(sprintf("Position: %s\n", format(pos, big.mark=",")))
      cat(sprintf("Sample: %s\n", sample_name))
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
        cat(sprintf("\n--- Window %d: %s ---\n", window_idx, 
                    ifelse(window_size >= 1000, paste0(window_size/1000, " kb"), paste0(window_size, " bp"))))
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
        cat(sprintf("SNPs in window: %d positions\n", nrow(founder_matrix_clean)))
        cat(sprintf("Clustering results: %d groups\n", n_groups))
        
        if (verbose >= 3 && n_groups <= 6) {
          # Show group assignments
          unique_clusters <- unique(groups)
          for (cluster_id in unique_clusters) {
            cluster_founders <- founders[groups == cluster_id]
            if (length(cluster_founders) == 1) {
              cat(sprintf("  Group %d: %s (individual)\n", cluster_id, cluster_founders))
            } else {
              cat(sprintf("  Group %d: %s\n", cluster_id, paste(cluster_founders, collapse=", ")))
            }
          }
        }
      }
      
      # Check if clustering improved
      if (window_idx > 1 && n_groups <= previous_n_groups) {
        if (verbose >= 2) cat("✗ No improvement from previous window, trying larger window\n")
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
            cat(sprintf("    Building constraints from %d founder groups...\n", n_groups))
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
              
              if (verbose >= 2) {
                group_names <- founders[cluster_founders]
                cat(sprintf("      Group %d: %s = %.4f\n", 
                           cluster_id, paste(group_names, collapse="+"), group_freq))
              }
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
              # Create constraint: this founder = their exact frequency
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
              
              if (verbose >= 2) {
                cat(sprintf("      Group %d: %s = %.4f (locked)\n", 
                           cluster_id, founders[cluster_founders], founder_freq))
              }
            }
          }
          
          # Update accumulated constraints for next window
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
            if (verbose >= 2) {
              cat(sprintf("    ✓ Accumulated %d constraints for next window (groups: %s)\n", 
                         nrow(current_constraints),
                         paste(unique(groups), collapse=",")))
            }
          } else {
            accumulated_constraints <- NULL
            accumulated_constraint_values <- NULL
            if (verbose >= 2) {
              cat("    • No constraints to accumulate (all founders in 1 group)\n")
            }
          }
          
          # Store the result (will be overwritten as we expand)
          final_result <- result
          final_n_groups <- n_groups
          
          if (verbose >= 2) {
            cat(sprintf("    Stored LSEI result with %d groups\n", n_groups))
            if ("cov" %in% names(result)) {
              cat("    Covariance matrix available in this result\n")
            } else {
              cat("    No covariance matrix in this result\n")
            }
          }
          
          # Check if all founders are separated
          if (n_groups == length(founders)) {
            if (verbose >= 2) {
              cat("✅ SUCCESS: All founders distinguished - stopping window expansion\n")
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
        if (verbose >= 2) {
          cat("✓ Founders distinguishable - frequencies are TRUSTWORTHY (estimate_OK = 1)\n")
        }
      } else {
        estimate_OK <- 0  # Founders NOT distinguishable - we CANNOT trust the frequencies
        if (verbose >= 2) {
          cat("⚠️  Founders NOT fully distinguished - frequencies are UNTRUSTWORTHY (estimate_OK = 0)\n")
        }
      }
      
    } else {
      # Either insufficient SNPs OR LSEI failed/didn't converge
      founder_frequencies <- rep(NA, length(founders))
      names(founder_frequencies) <- founders
      estimate_OK <- NA
      if (verbose >= 1) {
        cat("❌ LSEI failed or insufficient SNPs - no haplotype frequencies available\n")
      }
    }
    
    if (verbose >= 2) {
      cat(sprintf("\n✅ FINAL RESULT:\n"))
      cat(sprintf("  estimate_OK: %s\n", ifelse(is.na(estimate_OK), "NA", estimate_OK)))
      cat(sprintf("  final_window_size: %s\n", 
                  ifelse(final_window_size >= 1000, paste0(final_window_size/1000, " kb"), paste0(final_window_size, " bp"))))
      cat(sprintf("  n_snps: %d\n", ifelse(is.null(final_result), 0, nrow(final_wide_data))))
    }
    
    # Return appropriate error matrix
    if (!is.null(final_result)) {
      # Check what fields are available in the LSEI result
      if (verbose >= 2) {
        cat("LSEI result fields:", names(final_result), "\n")
        if ("cov" %in% names(final_result)) {
          cat("Covariance matrix available\n")
        } else {
          cat("No covariance matrix in LSEI result\n")
        }
      }
      
      if ("cov" %in% names(final_result) && !is.null(final_result$cov)) {
        error_matrix <- final_result$cov
        if (verbose >= 2) {
          cat("Using LSEI covariance matrix\n")
        }
      } else {
        # LSEI didn't provide covariance - create NA matrix
        error_matrix <- matrix(NA, length(founders), length(founders))
        rownames(error_matrix) <- founders
        colnames(error_matrix) <- founders
        if (verbose >= 2) {
          cat("LSEI result fields:", paste(names(final_result), collapse=", "), "\n")
          cat("No covariance matrix in LSEI result - using NA matrix\n")
        }
      }
    } else {
      # Create NA matrix when no result
      error_matrix <- matrix(NA, length(founders), length(founders))
      rownames(error_matrix) <- founders
      colnames(error_matrix) <- founders
      if (verbose >= 2) {
        cat("No LSEI result - using NA matrix\n")
      }
    }
    
    return(list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders))
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

run_adaptive_estimation <- function(chr, method, parameter, output_dir, param_file, debug = FALSE, verbose = TRUE) {
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
    all_positions <- head(all_positions, 100)  # Limit to 100 positions for debugging
    names_in_bam <- head(names_in_bam, 1)      # Limit to 1 sample for debugging
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
          verbose = ifelse(verbose, 1, 0)
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
          verbose = 0  # Always non-verbose for large datasets
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
    stop("Usage: Rscript complete_haplotype_workflow.R <chr> <method> <parameter> <output_dir> <param_file> [--nonverbose]")
  }
  
  chr <- args[1]
  method <- args[2]
  parameter <- as.numeric(args[3])
  output_dir <- args[4]
  param_file <- args[5]
  debug <- "--debug" %in% args
  verbose <- !("--nonverbose" %in% args)
  
  # Run the complete workflow
  if (debug) {
    cat("=== COMPLETE HAPLOTYPE WORKFLOW (DEBUG MODE) ===\n")
    cat("Limited to 100 positions × 1 sample for testing\n")
  } else {
    cat("=== COMPLETE HAPLOTYPE WORKFLOW ===\n")
    cat("Processing all positions and samples\n")
  }
  
  if (verbose) {
    cat("Verbose output enabled\n\n")
  } else {
    cat("Minimal output mode\n\n")
  }
  
  # Step 1: Adaptive estimation
  adaptive_results <- run_adaptive_estimation(chr, method, parameter, output_dir, param_file, debug, verbose)
  
  # Step 2: Smoothing
  smooth_results <- run_smoothing(chr, param_file, output_dir, adaptive_results, verbose)
  
  cat("\n=== WORKFLOW COMPLETE ===\n")
  cat("✓ Both adaptive estimation and smoothing completed successfully\n")
  cat("✓ All output files created in the same format as production\n")
}