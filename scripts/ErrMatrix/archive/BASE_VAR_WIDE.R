#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(limSolve)
library(MASS)

# =============================================================================
# BASE_VAR_WIDE.R - Wide Format Optimized Haplotype Estimation
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
#   PURPOSE: Core adaptive window haplotype estimation function
#   INPUT:  - pos: test position (bp)
#           - sample_name: sample to estimate haplotypes for
#           - df3: long-format data (POS, name, freq)
#           - founders: vector of founder names
#           - h_cutoff: clustering threshold (typically 4)
#   OUTPUT: List with Groups, Haps, Err, Names (HARDWIRED format)
#   LOGIC:  - Tests 6 window sizes: 10kb, 20kb, 50kb, 100kb, 200kb, 500kb
#           - For each window: filters data, pivots to wide, clusters founders
#           - Runs LSEI with accumulated constraints from previous windows
#           - Stops when all founders are distinguishable (8 groups)
#           - Returns frequency estimates and error matrix
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
# 3. run_smoothing() → takes adaptive results, applies 21-position smoothing
# 4. Main execution → runs both steps in sequence
#
# DATA FLOW:
# ----------
# RefAlt.txt → process_refalt_data() → df3 (long format)
# df3 → estimate_haplotypes_list_format() → Groups, Haps, Err, Names
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
# PROGRESSIVE ERROR MATRIX VERSION - SWAPPED IN
# =============================================================================

# Production function for delegation
estimate_haplotypes_list_format_prod <- function(pos, sample_name, df3, founders, h_cutoff,
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
  
  # No redundant debug output - function is called for every position
  
  method <- match.arg(method)
  
  # No redundant position processing output - already shown by wrapper
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
              
              # Constraint building details removed for cleaner output
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
              # Create constraint: this founder = their exact frequency
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
              
              # Constraint building details removed for cleaner output
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
          
          if (verbose >= 2) {
            cat(sprintf("    Stored LSEI result with %d groups\n", n_groups))
            if ("covar" %in% names(result)) {
              cat("    Covariance matrix available in this result (covar field)\n")
            } else {
              cat("    No covariance matrix in this result\n")
            }
          }
          
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
        if ("covar" %in% names(final_result)) {
          cat("Covariance matrix available (covar field)\n")
        } else {
          cat("No covariance matrix in LSEI result\n")
        }
      }
      
      if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
        error_matrix <- final_result$covar
        if (verbose >= 2) {
          cat("Using LSEI covariance matrix (covar field)\n")
        }
        if (verbose >= 2) {
          cat("\nError matrix (×100, rounded to 2 decimals):\n")
          print(round(error_matrix * 100, 2))
          cat("Condition number:", round(kappa(error_matrix), 2), "\n")
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

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {

  if (verbose >= 2) {
    cat(sprintf("=== PROGRESSIVE ERROR MATRIX MODE: h_cutoff=%g, sample=%s, pos=%g ===\n", 
                h_cutoff, sample_name, pos))
  }

  window_sizes <- c(150, 300, 750, 1500, 3000)
  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founders))
  groups_path <- integer(0)
  groups_fmt_path <- character(0)
  groups_win_path <- integer(0)

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  # Progressive V matrix (8x8), pending covariances, and resolved tracker
  V <- matrix(NA_real_, length(founders), length(founders))
  rownames(V) <- founders; colnames(V) <- founders
  pending <- list()
  resolved <- setNames(rep(FALSE, length(founders)), founders)

  # Track last (largest successful) unpooled design for true-cov computation
  last_A <- NULL
  last_y <- NULL

  # Helpers to compute pooled covariance at a window and to print V compactly
  pooled_cov <- function(A_full, y_full, groups_vec){
    gid <- sort(unique(groups_vec))
    k <- length(gid)
    A_pool <- matrix(0, nrow(A_full), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups_vec == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
    }
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar)) {
      # #tempdebug - Check what lsei returned
      cat("  lsei returned covar class:", class(fit$covar), "\n")
      cat("  lsei returned covar dimensions:", dim(fit$covar), "\n")
      if (!is.matrix(fit$covar)) {
        cat("  Converting lsei covar to matrix\n")
        fit$covar <- matrix(fit$covar, nrow=k, ncol=k)
      }
      return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
    }
    XtX <- crossprod(A_pool)
    xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) ginv(XtX) %*% crossprod(A_pool, y_full))
    r <- y_full - as.numeric(A_pool %*% xhat)
    sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
    cov_matrix <- sigma2 * tryCatch(solve(XtX), error=function(e) ginv(XtX))
    # Ensure cov is always a matrix, even for single group
    if (!is.matrix(cov_matrix)) {
      cov_matrix <- matrix(cov_matrix, nrow=k, ncol=k)
    }
    list(cov = cov_matrix, members=members, w=as.numeric(xhat))
  }

  fmt_cell_signed <- function(x, diag_cell){
    if (is.na(x)) return("  NA ")
    if (x == 0) {
      signc <- if (diag_cell) " " else "+"
      return(paste0(signc, "  0 "))
    }
    signc <- if (diag_cell) " " else if (x < 0) "-" else "+"
    ax <- abs(x)
    expo <- floor(log10(ax))
    mant <- ax / (10^expo)
    m2 <- as.integer(round(mant * 10))  # two-digit mantissa (approx)
    if (expo < -9) expo <- -9
    if (expo > 9) expo <- 9
    paste0(signc, sprintf("%02d%+1d", m2, as.integer(expo)))
  }

  print_V_compact <- function(Vmat){
    cat("V (signed covariances, sci 2d+exp; diag no sign):\n")
    for (i in seq_len(nrow(Vmat))){
      row <- character(ncol(Vmat))
      for (j in seq_len(ncol(Vmat))){
        row[j] <- fmt_cell_signed(Vmat[i, j], diag_cell = (i == j))
      }
      cat(paste(row, collapse=" "), "\n", sep="")
    }
  }

  # Helpers to compute pooled covariance at a window and to print V compactly
  pooled_cov <- function(A_full, y_full, groups_vec){
    gid <- sort(unique(groups_vec))
    k <- length(gid)
    A_pool <- matrix(0, nrow(A_full), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups_vec == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
    }
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar)) {
      # #tempdebug - Check what lsei returned
      cat("  lsei returned covar class:", class(fit$covar), "\n")
      cat("  lsei returned covar dimensions:", dim(fit$covar), "\n")
      if (!is.matrix(fit$covar)) {
        cat("  Converting lsei covar to matrix\n")
        fit$covar <- matrix(fit$covar, nrow=k, ncol=k)
      }
      return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
    }
    XtX <- crossprod(A_pool)
    xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) ginv(XtX) %*% crossprod(A_pool, y_full))
    r <- y_full - as.numeric(A_pool %*% xhat)
    sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
    cov_matrix <- sigma2 * tryCatch(solve(XtX), error=function(e) ginv(XtX))
    # Ensure cov is always a matrix, even for single group
    if (!is.matrix(cov_matrix)) {
      cov_matrix <- matrix(cov_matrix, nrow=k, ncol=k)
    }
    list(cov = cov_matrix, members=members, w=as.numeric(xhat))
  }

  fmt_cell_signed <- function(x, diag_cell){
    if (is.na(x)) return("  NA ")
    if (x == 0) {
      signc <- if (diag_cell) " " else "+"
      return(paste0(signc, "  0 "))
    }
    signc <- if (diag_cell) " " else if (x < 0) "-" else "+"
    ax <- abs(x)
    expo <- floor(log10(ax))
    mant <- ax / (10^expo)
    m2 <- as.integer(round(mant * 10))  # two-digit mantissa (approx)
    if (expo < -9) expo <- -9
    if (expo > 9) expo <- 9
    paste0(signc, sprintf("%02d%+1d", m2, as.integer(expo)))
  }

  print_V_compact <- function(Vmat){
    cat("V (signed covariances, sci 2d+exp; diag no sign):\n")
    for (i in seq_len(nrow(Vmat))){
      row <- character(ncol(Vmat))
      for (j in seq_len(ncol(Vmat))){
        row[j] <- fmt_cell_signed(Vmat[i, j], diag_cell = (i == j))
      }
      cat(paste(row, collapse=" "), "\n", sep="")
    }
  }

  # Ensure POS-ordered tidy input (POS, name, freq)
  df3 <- df3 %>% dplyr::arrange(POS)
  
  for (window_size in window_sizes) {
    # Convert SNP count to base pairs (approximate: 1 SNP per 100bp)
    window_size_bp <- window_size * 100
    window_start <- max(1, pos - window_size_bp/2)
    window_end <- pos + window_size_bp/2
    
    window_data <- df3 %>%
      dplyr::filter(POS >= window_start, POS <= window_end, name %in% c(founders, sample_name))
    if (nrow(window_data) == 0) next

    wide_data <- window_data %>%
      dplyr::select(POS, name, freq) %>%
      tidyr::pivot_wider(names_from = name, values_from = freq)

    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next

    founder_matrix <- wide_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
    sample_freqs <- wide_data %>% dplyr::pull(!!sample_name)

    complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    if (nrow(founder_matrix_clean) < 10) next

    # Clustering (Euclidean, complete linkage), cutree at h_cutoff
    founder_dist <- stats::dist(t(founder_matrix_clean))
    hclust_result <- stats::hclust(founder_dist, method = "complete")
    groups <- stats::cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) {
      comp <- split(founders, groups)
      ordered_gids <- sort(as.integer(names(comp)))
      group_strings <- vapply(ordered_gids, function(gid){
        paste(comp[[as.character(gid)]], collapse="+")
      }, character(1))
      carried_ct <- if (!is.null(accumulated_constraints)) nrow(accumulated_constraints) else 0
    }

    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("    No improvement; continue\n")
      next
    }
    previous_n_groups <- n_groups

    # Constraints accumulate like production
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    Carried over %d constraints\n", nrow(accumulated_constraints)))
    }

    res <- tryCatch(
      limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean,
                     E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("    LSEI failed; continue\n")
      next
    }

    # Build/accumulate group constraints (pool groups, lock singles)
    current_constraints <- NULL
    current_constraint_values <- NULL
    for (cid in unique(groups)) {
      idx <- which(groups == cid)
      row <- rep(0, n_founders); row[idx] <- 1
      current_constraints <- rbind(current_constraints, row)
      current_constraint_values <- c(current_constraint_values, sum(res$X[idx]))
    }
    if (!is.null(current_constraints)) {
      accumulated_constraints <- current_constraints
      accumulated_constraint_values <- current_constraint_values
      if (verbose >= 2) {
        built_ct <- nrow(current_constraints)
        # Map values to ordered group ids
        ordered_gids <- sort(unique(groups))
        group_values <- vapply(ordered_gids, function(gid){
          sum(res$X[which(groups == gid)])
        }, numeric(1))
      }
    } else {
      accumulated_constraints <- NULL
      accumulated_constraint_values <- NULL
      if (verbose >= 2) {
        built_ct <- 0
        group_values <- numeric(0)
      }
    }

    # Update progressive V using pooled model for this window
    pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
    cov_pool <- pc$cov; pool_members <- pc$members
    
    # Debug: check what we got
    if (verbose >= 1) {
      cat("cov_pool class:", class(cov_pool), "\n")
      cat("cov_pool dimensions:", dim(cov_pool), "\n")
      cat("pool_members length:", length(pool_members), "\n")
      cat("groups:", groups, "\n")
    }
    
    # Safety check: ensure cov_pool is a matrix
    if (!is.matrix(cov_pool) || any(is.na(cov_pool))) {
      if (verbose >= 1) {
        cat("Warning: cov_pool is not a valid matrix, skipping progressive V update\n")
      }
      next
    }
    # mark newly resolved founders
    for (ii in seq_along(pool_members)){
      mem <- pool_members[[ii]]
      if (length(mem)==1L) resolved[founders[mem]] <- TRUE
    }
    # write cov for resolved-resolved; store pending for det-vs-pool and pool-vs-pool
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        if (length(mi)==1L && length(mj)==1L){
          fi <- founders[mi]; fj <- founders[mj]
          if (is.na(V[fi,fj]) && i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            V[fi,fj] <- cov_pool[i_idx, j_idx]
          }
          if (is.na(V[fj,fi]) && j_idx <= nrow(cov_pool) && i_idx <= ncol(cov_pool)) {
            V[fj,fi] <- cov_pool[j_idx, i_idx]
          }
        } else if (length(mi)==1L && length(mj)>1L){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="det_vs_pool", det=founders[mi], pool=j_idx, cov=cov_pool[i_idx, j_idx], mem=founders[mj])))
          }
        } else if (length(mi)>1L && length(mj)>1L && i_idx<=j_idx){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="pool_vs_pool", pool_a=i_idx, pool_b=j_idx, cov=cov_pool[i_idx, j_idx], mem_a=founders[mi], mem_b=founders[mj])))
          }
        }
      }
    }

    # If all resolved now, distribute pending covariances using current weights
    if (all(resolved)){
      w_now <- setNames(rep(NA_real_, length(founders)), founders)
      # Approximate weights from current (last) res if same dimensionality, else use group_values map
      if (length(final_result$X) == length(founders)){
        w_now <- setNames(as.numeric(final_result$X), founders)
      }
      for (rec in pending){
        if (rec$type == "det_vs_pool"){
          den <- sum(w_now[rec$mem]); if (!is.finite(den) || den<=0) next
          for (f in rec$mem){
            val <- rec$cov * (w_now[f]/den)
            if (is.na(V[rec$det, f])) V[rec$det, f] <- val
            if (is.na(V[f, rec$det])) V[f, rec$det] <- val
          }
        } else if (rec$type == "pool_vs_pool"){
          den_a <- sum(w_now[rec$mem_a]); den_b <- sum(w_now[rec$mem_b]); if (den_a<=0 || den_b<=0) next
          for (fa in rec$mem_a){
            for (fb in rec$mem_b){
              val <- rec$cov * (w_now[fa]/den_a) * (w_now[fb]/den_b)
              if (is.na(V[fa, fb])) V[fa, fb] <- val
              if (is.na(V[fb, fa])) V[fb, fa] <- val
            }
          }
        }
      }
      pending <- list()
    }

    if (verbose >= 2){
      print_V_compact(V)
    }
    
    if (n_groups == length(founders)) {
      # Save unpooled design at the final successful window
      last_A <- founder_matrix_clean
      last_y <- sample_freqs_clean
      if (verbose >= 2) cat("    All founders distinguished; stop\n")
      break
    }

    final_result <- res
    final_n_groups <- n_groups
    groups_path <- c(groups_path, n_groups)
    # compact groups format like [123]4[56]78 based on founders order
    fmt <- {
      grp_ids <- sort(unique(groups))
      parts <- character(0)
      for (gid in grp_ids) {
        idx <- which(groups == gid)
        if (length(idx) > 1) {
          parts <- c(parts, paste0("[", paste(idx, collapse=""), "]"))
        } else {
          parts <- c(parts, paste0(idx))
        }
      }
      paste(parts, collapse = "")
    }
    groups_fmt_path <- c(groups_fmt_path, fmt)
    groups_win_path <- c(groups_win_path, window_size)
    if (verbose >= 2) {
      # Print compact, aligned 3-line block (monospace-friendly)
      cat(sprintf("window_snp=%-5d  n_snps=%-5d  n_groups=%-2d  carried=%-2d  built=%-2d\n",
                  window_size, nrow(founder_matrix_clean), n_groups, carried_ct, built_ct))
      if (length(group_strings)) {
        # Column widths match label lengths exactly for alignment
        col_w <- nchar(group_strings)
        # Render group labels with padding
        grp_fmt <- vapply(seq_along(group_strings), function(i){
          sprintf("%-*s", col_w[i], group_strings[i])
        }, character(1))
        cat(paste(grp_fmt, collapse=" | "), "\n", sep="")
      }
      if (length(group_values)) {
        # Convert to percent, 0 decimals
        vals_pct <- round(group_values * 100, 0)
        # Match widths of labels exactly
        if (!exists("col_w")) col_w <- rep(1, length(vals_pct))
        val_fmt <- vapply(seq_along(vals_pct), function(i){
          sprintf("%*d", col_w[i], as.integer(vals_pct[i]))
        }, character(1))
        cat(paste(val_fmt, collapse=" | "), "\n", sep="")
      }
    }

    # Update progressive V using pooled model for this window
    pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
    cov_pool <- pc$cov; pool_members <- pc$members
    
    # Debug: check what we got
    if (verbose >= 1) {
      cat("cov_pool class:", class(cov_pool), "\n")
      cat("cov_pool dimensions:", dim(cov_pool), "\n")
      cat("pool_members length:", length(pool_members), "\n")
      cat("groups:", groups, "\n")
    }
    
    # Safety check: ensure cov_pool is a matrix
    if (!is.matrix(cov_pool) || any(is.na(cov_pool))) {
      if (verbose >= 1) {
        cat("Warning: cov_pool is not a valid matrix, skipping progressive V update\n")
      }
      next
    }
    # mark newly resolved founders
    for (ii in seq_along(pool_members)){
      mem <- pool_members[[ii]]
      if (length(mem)==1L) resolved[founders[mem]] <- TRUE
    }
    # write cov for resolved-resolved; store pending for det-vs-pool and pool-vs-pool
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        if (length(mi)==1L && length(mj)==1L){
          fi <- founders[mi]; fj <- founders[mj]
          if (is.na(V[fi,fj]) && i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            V[fi,fj] <- cov_pool[i_idx, j_idx]
          }
          if (is.na(V[fj,fi]) && j_idx <= nrow(cov_pool) && i_idx <= ncol(cov_pool)) {
            V[fj,fi] <- cov_pool[j_idx, i_idx]
          }
        } else if (length(mi)==1L && length(mj)>1L){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="det_vs_pool", det=founders[mi], pool=j_idx, cov=cov_pool[i_idx, j_idx], mem=founders[mj])))
          }
        } else if (length(mi)>1L && length(mj)>1L && i_idx<=j_idx){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="pool_vs_pool", pool_a=i_idx, pool_b=j_idx, cov=cov_pool[i_idx, j_idx], mem_a=founders[mi], mem_b=founders[mj])))
          }
        }
      }
    }

    # If all resolved now, distribute pending covariances using current weights
    if (all(resolved)){
      w_now <- setNames(rep(NA_real_, length(founders)), founders)
      # Approximate weights from current (last) res if same dimensionality, else use group_values map
      if (length(final_result$X) == length(founders)){
        w_now <- setNames(as.numeric(final_result$X), founders)
      }
      for (rec in pending){
        if (rec$type == "det_vs_pool"){
          den <- sum(w_now[rec$mem]); if (!is.finite(den) || den<=0) next
          for (f in rec$mem){
            val <- rec$cov * (w_now[f]/den)
            if (is.na(V[rec$det, f])) V[rec$det, f] <- val
            if (is.na(V[f, rec$det])) V[f, rec$det] <- val
          }
        } else if (rec$type == "pool_vs_pool"){
          den_a <- sum(w_now[rec$mem_a]); den_b <- sum(w_now[rec$mem_b]); if (den_a<=0 || den_b<=0) next
          for (fa in rec$mem_a){
            for (fb in rec$mem_b){
              val <- rec$cov * (w_now[fa]/den_a) * (w_now[fb]/den_b)
              if (is.na(V[fa, fb])) V[fa, fb] <- val
              if (is.na(V[fb, fa])) V[fb, fa] <- val
            }
          }
        }
      }
      pending <- list()
    }


    if (verbose >= 2) {
      cat("V matrix:\n")
      print_V_compact(V)
    }

    # Check if all founders are resolved
    if (all(resolved)) {
      if (verbose >= 2) cat("    All founders resolved - stopping\n")
      break
    }
  }

  # Return results
  if (is.null(final_result)) {
    return(NULL)
  }
  
  # Calculate founder frequencies and error matrix
  founder_frequencies <- final_result$X
  names(founder_frequencies) <- founders
  
  # Compute "true" covariance using largest window residual sigma^2 and unpooled design
  Cov_true <- NULL
  if (!is.null(last_A) && !is.null(last_y)) {
    XtX <- crossprod(last_A)
    Xinv <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
    r <- last_y - as.numeric(last_A %*% founder_frequencies)
    p <- ncol(last_A); n <- nrow(last_A)
    sigma2_hat <- sum(r^2) / max(1, n - p)
    Cov_true <- sigma2_hat * Xinv
  }
  
  # Prefer progressive V if sufficiently filled; otherwise fallback
  if (sum(is.na(V)) < length(V)) {
    for (d in seq_len(nrow(V))) if (is.na(V[d, d])) V[d, d] <- 1e-8
    error_matrix <- V
  } else if (!is.null(Cov_true)) {
    error_matrix <- Cov_true
  } else if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
    error_matrix <- final_result$covar
  } else {
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
  }
  
  res_out <- list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders)
  # Attach diagnostics for wrapper comparison
  attr(res_out, "V_progressive") <- V
  if (exists("Cov_true")) attr(res_out, "Cov_true") <- Cov_true else attr(res_out, "Cov_true") <- NULL
  attr(res_out, "groups_path") <- groups_path
  attr(res_out, "groups_fmt_path") <- groups_fmt_path
  attr(res_out, "groups_win_path") <- groups_win_path
  res_out
}

# =============================================================================
# DATA PROCESSING FUNCTIONS (EXACT FROM 49H AGO)
# =============================================================================

process_refalt_data_wide <- function(refalt_file, founders) {
  cat("Loading RefAlt data in WIDE format from:", refalt_file, "\n")
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  cat("Raw RefAlt data rows:", nrow(refalt_data), "\n")
  
  # Keep data in wide format - no pivot_longer!
  # Extract sample names and create frequency matrices
  sample_names <- unique(str_sub(names(refalt_data)[-c(1,2)], 5))  # Remove CHROM, POS, extract sample names
  
  # Create REF and ALT matrices for each sample
  ref_cols <- names(refalt_data)[str_sub(names(refalt_data), 1, 3) == "REF" & str_sub(names(refalt_data), 5) %in% sample_names]
  alt_cols <- names(refalt_data)[str_sub(names(refalt_data), 1, 3) == "ALT" & str_sub(names(refalt_data), 5) %in% sample_names]
  
  ref_matrix <- as.matrix(refalt_data[, ref_cols])
  alt_matrix <- as.matrix(refalt_data[, alt_cols])
  
  # Calculate frequencies and read counts
  freq_matrix <- ref_matrix / (ref_matrix + alt_matrix)
  n_matrix <- ref_matrix + alt_matrix
  
  # Set column names to sample names
  colnames(freq_matrix) <- sample_names
  colnames(n_matrix) <- sample_names
  
  # Quality filtering: EXACT same as BASE_VAR.R - keep positions where ALL founders are "fixed"
  founder_freqs <- freq_matrix[, founders, drop = FALSE]
  valid_positions <- apply(founder_freqs, 1, function(row) {
    all(is.na(row) | row < 0.03 | row > 0.97)
  })
  
  # Filter to valid positions
  cat("Valid positions after filtering:", sum(valid_positions), "out of", length(valid_positions), "\n")
  refalt_wide <- refalt_data[valid_positions, ]
  freq_matrix <- freq_matrix[valid_positions, ]
  n_matrix <- n_matrix[valid_positions, ]
  
  # Set column names
  colnames(freq_matrix) <- sample_names
  colnames(n_matrix) <- sample_names
  
  cat("✓ Processed", nrow(refalt_wide), "positions for", length(sample_names), "samples in WIDE format\n")
  
  return(list(
    refalt = refalt_wide,
    freq_matrix = freq_matrix,
    n_matrix = n_matrix,
    sample_names = sample_names
  ))
}

# Wide format haplotype estimation function
estimate_haplotypes_wide <- function(pos, sample_name, wide_data, founders, h_cutoff,
                                   method = "adaptive",
                                   window_size_bp = NULL,
                                   chr = "chr2R",
                                   verbose = 0) {
  
  if (verbose >= 2) {
    cat(sprintf("=== WIDE FORMAT HAPLOTYPE ESTIMATION: h_cutoff=%g, sample=%s, pos=%g ===\n", 
                h_cutoff, sample_name, pos))
  }

  # Work directly with wide format matrices - this is the key optimization!
  # We avoid the expensive pivot_longer/pivot_wider operations entirely
  positions <- wide_data$refalt$POS
  freq_matrix <- wide_data$freq_matrix
  
  # Call the working haplotype estimation function but pass the wide data
  # The function will need to be modified to accept wide format data
  return(estimate_haplotypes_list_format_wide_input(pos, sample_name, freq_matrix, positions, founders, h_cutoff, method, window_size_bp, chr, verbose))
}

# Wide format input version of the haplotype estimation function
estimate_haplotypes_list_format_wide_input <- function(pos, sample_name, freq_matrix, positions, founders, h_cutoff,
                                   method = "adaptive",
                                   window_size_bp = NULL,
                                   chr = "chr2R",
                                   verbose = 0) {

  if (verbose >= 2) {
    cat(sprintf("=== PROGRESSIVE ERROR MATRIX MODE: h_cutoff=%g, sample=%s, pos=%g ===\n",
                h_cutoff, sample_name, pos))
  }

  window_sizes <- c(150, 300, 750, 1500, 3000) # These are SNP counts, converted to BP
  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founders))
  groups_path <- integer(0)
  groups_fmt_path <- character(0)
  groups_win_path <- integer(0)

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  # Progressive V matrix (8x8), pending covariances, and resolved tracker
  V <- matrix(NA_real_, length(founders), length(founders))
  rownames(V) <- founders; colnames(V) <- founders
  pending <- list()
  resolved <- setNames(rep(FALSE, length(founders)), founders)

  # Track last (largest successful) unpooled design for true-cov computation
  last_A <- NULL
  last_y <- NULL

  # Helpers to compute pooled covariance at a window and to print V compactly
  pooled_cov <- function(A_full, y_full, groups_vec){
    gid <- sort(unique(groups_vec))
    k <- length(gid)
    A_pool <- matrix(0, nrow(A_full), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups_vec == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
    }
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar)) {
      # #tempdebug - Check what lsei returned
      cat("  lsei returned covar class:", class(fit$covar), "\n")
      cat("  lsei returned covar dimensions:", dim(fit$covar), "\n")
      if (!is.matrix(fit$covar)) {
        cat("  Converting lsei covar to matrix\n")
        fit$covar <- matrix(fit$covar, nrow=k, ncol=k)
      }
      return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
    }
    XtX <- crossprod(A_pool)
    xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) ginv(XtX) %*% crossprod(A_pool, y_full))
    r <- y_full - as.numeric(A_pool %*% xhat)
    sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
    cov_matrix <- sigma2 * tryCatch(solve(XtX), error=function(e) ginv(XtX))
    # Ensure cov is always a matrix, even for single group
    if (!is.matrix(cov_matrix)) {
      cov_matrix <- matrix(cov_matrix, nrow=k, ncol=k)
    }
    list(cov = cov_matrix, members=members, w=as.numeric(xhat))
  }

  fmt_cell_signed <- function(x, diag_cell){
    if (is.na(x)) return("  NA ")
    if (x == 0) {
      signc <- if (diag_cell) " " else "+"
      return(paste0(signc, "  0 "))
    }
    signc <- if (diag_cell) " " else if (x < 0) "-" else "+"
    ax <- abs(x)
    expo <- floor(log10(ax))
    mant <- ax / (10^expo)
    m2 <- as.integer(round(mant * 10))
    if (expo < -9) expo <- -9
    if (expo > 9) expo <- 9
    paste0(signc, sprintf("%02d%+1d", m2, as.integer(expo)))
  }

  print_V_compact <- function(Vmat){
    cat("V (signed covariances, sci 2d+exp; diag no sign):\n")
    for (i in seq_len(nrow(Vmat))){
      row <- character(ncol(Vmat))
      for (j in seq_len(ncol(Vmat))){
        row[j] <- fmt_cell_signed(Vmat[i, j], diag_cell = (i == j))
      }
      cat(paste(row, collapse=" "), "\n", sep="")
    }
  }

  # #tempdebug - Print input data to reality check
  cat("=== INPUT DATA DEBUG ===\n")
  cat("freq_matrix class:", class(freq_matrix), "\n")
  cat("freq_matrix dimensions:", dim(freq_matrix), "\n")
  cat("freq_matrix column names:", paste(colnames(freq_matrix), collapse=", "), "\n")
  cat("positions length:", length(positions), "\n")
  cat("positions range:", range(positions), "\n")
  cat("sample_name:", sample_name, "\n")
  cat("founders:", paste(founders, collapse=", "), "\n")
  cat("First few rows of freq_matrix:\n")
  print(head(freq_matrix, 3))
  cat("========================\n")
  
  # #tempdebug - Count positions for first estimation
  if (pos == 5400000 && sample_name == "AJ_1_1") {
    cat("=== BASE_VAR_WIDE DEBUG: First estimation ===\n")
    cat("Total freq_matrix rows:", nrow(freq_matrix), "\n")
    cat("freq_matrix columns:", paste(colnames(freq_matrix), collapse=", "), "\n")
    cat("positions range:", range(positions), "\n")
    cat("Sample:", sample_name, "Founders:", paste(founders, collapse=", "), "\n")
  }

  # Work directly with wide format matrices - no pivot operations needed!
  for (window_size_snp_count in window_sizes) {
    # Convert SNP count to base pairs (approximate: 1 SNP per 100bp)
    window_size_bp_actual <- window_size_snp_count * 100
    window_start <- max(1, pos - window_size_bp_actual/2)
    window_end <- pos + window_size_bp_actual/2

    # Find positions in window
    window_positions <- positions >= window_start & positions <= window_end
    if (sum(window_positions) == 0) next
    
    # Extract window data directly from matrices
    window_freq_matrix <- freq_matrix[window_positions, , drop = FALSE]
    
    # #tempdebug - Print window data
    cat(sprintf("Window %d SNPs: %d positions, window_start=%d, window_end=%d\n", 
                window_size_snp_count, nrow(window_freq_matrix), window_start, window_end))
    cat("Window freq_matrix dimensions:", dim(window_freq_matrix), "\n")
    cat("Window freq_matrix column names:", paste(colnames(window_freq_matrix), collapse=", "), "\n")
    
    # #tempdebug - Count positions for first estimation, first window
    if (pos == 5400000 && sample_name == "AJ_1_1" && window_size_snp_count == 150) {
      cat("=== BASE_VAR_WIDE WINDOW DEBUG ===\n")
      cat("Window size:", window_size_snp_count, "SNPs\n")
      cat("Window start:", window_start, "end:", window_end, "\n")
      cat("Window freq_matrix rows:", nrow(window_freq_matrix), "\n")
      cat("Window positions range:", range(positions[window_positions]), "\n")
      cat("Window freq_matrix column names:", paste(colnames(window_freq_matrix), collapse=", "), "\n")
    }
    
    if (nrow(window_freq_matrix) < 10) {
      cat("  -> Skipping: too few positions\n")
      next
    }
    if (!all(c(founders, sample_name) %in% colnames(window_freq_matrix))) {
      cat("  -> Skipping: missing founders or sample\n")
      next
    }

    # Extract founder matrix and sample frequencies directly
    founder_matrix <- window_freq_matrix[, founders, drop = FALSE]
    sample_freqs <- window_freq_matrix[, sample_name]
    
    # #tempdebug - Print extracted data
    cat("  Founder matrix dimensions:", dim(founder_matrix), "\n")
    cat("  Sample freqs length:", length(sample_freqs), "\n")
    cat("  First few founder freqs:\n")
    print(head(founder_matrix, 3))
    cat("  First few sample freqs:", head(sample_freqs, 5), "\n")

    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    if (nrow(founder_matrix_clean) < 10) next

    # Clustering (Euclidean, complete linkage), cutree at h_cutoff
    founder_dist <- dist(t(founder_matrix_clean))
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) {
      comp <- split(founders, groups)
      ordered_gids <- sort(as.integer(names(comp)))
      group_strings <- vapply(ordered_gids, function(gid){
        paste(comp[[as.character(gid)]], collapse="+")
      }, character(1))
      carried_ct <- if (!is.null(accumulated_constraints)) nrow(accumulated_constraints) else 0
    }

    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("    No improvement; continue\n")
      next
    }
    previous_n_groups <- n_groups

    # Constraints accumulate like production
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    Carried over %d constraints\n", nrow(accumulated_constraints)))
    }

    res <- tryCatch(
      limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean,
                     E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("    LSEI failed; continue\n")
      next
    }

    # Build/accumulate group constraints (pool groups, lock singles)
    current_constraints <- NULL
    current_constraint_values <- NULL
    for (cid in unique(groups)) {
      idx <- which(groups == cid)
      row <- rep(0, n_founders); row[idx] <- 1
      current_constraints <- rbind(current_constraints, row)
      current_constraint_values <- c(current_constraint_values, sum(res$X[idx]))
    }
    if (!is.null(current_constraints)) {
      accumulated_constraints <- current_constraints
      accumulated_constraint_values <- current_constraint_values
      if (verbose >= 2) {
        built_ct <- nrow(current_constraints)
        # Map values to ordered group ids
        ordered_gids <- sort(unique(groups))
        group_values <- vapply(ordered_gids, function(gid){
          sum(res$X[which(groups == gid)])
        }, numeric(1))
      }
    } else {
      accumulated_constraints <- NULL
      accumulated_constraint_values <- NULL
      if (verbose >= 2) {
        built_ct <- 0
        group_values <- numeric(0)
      }
    }

    # Update progressive V using pooled model for this window
    pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
    cov_pool <- pc$cov; pool_members <- pc$members

    # Debug: check what we got
    if (verbose >= 1) {
      cat("cov_pool class:", class(cov_pool), "\n")
      cat("cov_pool dimensions:", dim(cov_pool), "\n")
      cat("pool_members length:", length(pool_members), "\n")
      cat("groups:", groups, "\n")
    }

    # Safety check: ensure cov_pool is a matrix
    if (!is.matrix(cov_pool) || any(is.na(cov_pool))) {
      if (verbose >= 1) {
        cat("Warning: cov_pool is not a valid matrix, skipping progressive V update\n")
      }
      next
    }

    # mark newly resolved founders
    for (ii in seq_along(pool_members)){
      mem <- pool_members[[ii]]
      if (length(mem)==1L) resolved[founders[mem]] <- TRUE
    }
    # write cov for resolved-resolved; store pending for det-vs-pool and pool-vs-pool
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        if (length(mi)==1L && length(mj)==1L){
          fi <- founders[mi]; fj <- founders[mj]
          if (is.na(V[fi,fj]) && i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            V[fi,fj] <- cov_pool[i_idx, j_idx]
          }
          if (is.na(V[fj,fi]) && j_idx <= nrow(cov_pool) && i_idx <= ncol(cov_pool)) {
            V[fj,fi] <- cov_pool[j_idx, i_idx]
          }
        } else if (length(mi)==1L && length(mj)>1L){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="det_vs_pool", det=founders[mi], pool=j_idx, cov=cov_pool[i_idx, j_idx], mem=founders[mj])))
          }
        } else if (length(mi)>1L && length(mj)>1L && i_idx<=j_idx){
          if (i_idx <= nrow(cov_pool) && j_idx <= ncol(cov_pool)) {
            pending <- append(pending, list(list(type="pool_vs_pool", pool_a=i_idx, pool_b=j_idx, cov=cov_pool[i_idx, j_idx], mem_a=founders[mi], mem_b=founders[mj])))
          }
        }
      }
    }

    # If all resolved now, distribute pending covariances using current weights
    if (all(resolved)){
      w_now <- setNames(rep(NA_real_, length(founders)), founders)
      # Approximate weights from current (last) res if same dimensionality, else use group_values map
      if (length(final_result$X) == length(founders)){
        w_now <- setNames(as.numeric(final_result$X), founders)
      }
      for (rec in pending){
        if (rec$type == "det_vs_pool"){
          den <- sum(w_now[rec$mem]); if (!is.finite(den) || den<=0) next
          for (f in rec$mem){
            val <- rec$cov * (w_now[f]/den)
            if (is.na(V[rec$det, f])) V[rec$det, f] <- val
            if (is.na(V[f, rec$det])) V[f, rec$det] <- val
          }
        } else if (rec$type == "pool_vs_pool"){
          den_a <- sum(w_now[rec$mem_a]); den_b <- sum(w_now[rec$mem_b]); if (den_a<=0 || den_b<=0) next
          for (fa in rec$mem_a){
            for (fb in rec$mem_b){
              val <- rec$cov * (w_now[fa]/den_a) * (w_now[fb]/den_b)
              if (is.na(V[fa, fb])) V[fa, fb] <- val
              if (is.na(V[fb, fa])) V[fb, fa] <- val
            }
          }
        }
      }
      pending <- list()
    }

    if (verbose >= 2) {
      cat("V matrix:\n")
      print_V_compact(V)
    }

    if (n_groups == length(founders)) {
      # Save unpooled design at the final successful window
      last_A <- founder_matrix_clean
      last_y <- sample_freqs_clean
      if (verbose >= 2) cat("    All founders distinguished; stop\n")
      break
    }

    final_result <- res
    final_n_groups <- n_groups
    groups_path <- c(groups_path, n_groups)
    # compact groups format like [123]4[56]78 based on founders order
    fmt <- {
      grp_ids <- sort(unique(groups))
      parts <- character(0)
      for (gid in grp_ids) {
        idx <- which(groups == gid)
        if (length(idx) > 1) {
          parts <- c(parts, paste0("[", paste(idx, collapse=""), "]"))
        } else {
          parts <- c(parts, paste0(idx))
        }
      }
      paste(parts, collapse = "")
    }
    groups_fmt_path <- c(groups_fmt_path, fmt)
    groups_win_path <- c(groups_win_path, window_size_snp_count) # Use original window_size for path
    if (verbose >= 2) {
      # Print compact, aligned 3-line block (monospace-friendly)
      cat(sprintf("window_snp=%-5d  n_snps=%-5d  n_groups=%-2d  carried=%-2d  built=%-2d\n",
                  window_size_snp_count, nrow(founder_matrix_clean), n_groups, carried_ct, built_ct))
      if (length(group_strings)) {
        # Column widths match label lengths exactly for alignment
        col_w <- nchar(group_strings)
        # Render group labels with padding
        grp_fmt <- vapply(seq_along(group_strings), function(i){
          sprintf("%-*s", col_w[i], group_strings[i])
        }, character(1))
        cat(paste(grp_fmt, collapse=" | "), "\n", sep="")
      }
      if (length(group_values)) {
        # Convert to percent, 0 decimals
        vals_pct <- round(group_values * 100, 0)
        # Match widths of labels exactly
        if (!exists("col_w")) col_w <- rep(1, length(vals_pct))
        val_fmt <- vapply(seq_along(vals_pct), function(i){
          sprintf("%*d", col_w[i], as.integer(vals_pct[i]))
        }, character(1))
        cat(paste(val_fmt, collapse=" | "), "\n", sep="")
      }
    }
  }

  # Return results
  if (is.null(final_result)) {
    return(NULL)
  }

  # Calculate founder frequencies and error matrix
  founder_frequencies <- final_result$X
  names(founder_frequencies) <- founders

  # Compute "true" covariance using largest window residual sigma^2 and unpooled design
  Cov_true <- NULL
  if (!is.null(last_A) && !is.null(last_y)) {
    XtX <- crossprod(last_A)
    Xinv <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
    r <- last_y - as.numeric(last_A %*% founder_frequencies)
    p <- ncol(last_A); n <- nrow(last_A)
    sigma2_hat <- sum(r^2) / max(1, n - p)
    Cov_true <- sigma2_hat * Xinv
  }

  # Prefer progressive V if sufficiently filled; otherwise fallback
  if (sum(is.na(V)) < length(V)) {
    for (d in seq_len(nrow(V))) if (is.na(V[d, d])) V[d, d] <- 1e-8
    error_matrix <- V
  } else if (!is.null(Cov_true)) {
    error_matrix <- Cov_true
  } else if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
    error_matrix <- final_result$covar
  } else {
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
  }

  # Check for singular matrix
  kappa_val <- kappa(error_matrix)
  if (is.infinite(kappa_val) || kappa_val > 1e12) {
    if (verbose >= 1) {
      cat("Warning: Singular error matrix (kappa =", kappa_val, "), using pseudoinverse\n")
    }
    error_matrix <- ginv(error_matrix)
  }

  res_out <- list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders)
  # Attach diagnostics for wrapper comparison
  attr(res_out, "V_progressive") <- V
  if (exists("Cov_true")) attr(res_out, "Cov_true") <- Cov_true else attr(res_out, "Cov_true") <- NULL
  attr(res_out, "groups_path") <- groups_path
  attr(res_out, "groups_fmt_path") <- groups_fmt_path
  attr(res_out, "groups_win_path") <- groups_win_path
  res_out
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
# WIDE FORMAT OPTIMIZED VERSION
# =============================================================================

run_adaptive_estimation_wide <- function(chr, method, parameter, output_dir, param_file, debug = FALSE, verbose = TRUE, debug_level = 0) {
  # Wide format optimized version - implements the optimization plan
  # 1. Read data in wide format (no pivot_longer)
  # 2. Step through estimation positions
  # 3. Subset the big dataframe for each position ± 250kb
  # 4. Map the estimation function to all samples for that window
  # 5. Pass wide data directly to the haplotype estimation function
  
  cat("=== RUNNING WIDE FORMAT ADAPTIVE HAPLOTYPE ESTIMATION ===\n")
  cat("Chromosome:", chr, "\n")
  cat("Method:", method, "\n")
  cat("Parameter:", parameter, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameter file:", param_file, "\n")
  cat("Debug mode:", debug, "\n")
  cat("Verbose mode:", verbose, "\n\n")
  
  # Load parameters
  source(param_file)
  
  # Load RefAlt data in WIDE format
  refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(refalt_file)) {
    stop("RefAlt file not found: ", refalt_file)
  }
  
  wide_data <- process_refalt_data_wide(refalt_file, founders)
  
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
  
  cat("Total positions to process:", length(all_positions), "\n")
  cat("Samples to process:", length(names_in_bam), "\n\n")
  
  # Initialize results storage
  results_list <- list()
  result_counter <- 0
  
  # Get the maximum window size (3000 SNPs = 300kb)
  max_window_bp <- 3000 * 100  # 300kb
  
  # Process each position
  for (pos_idx in seq_along(all_positions)) {
    pos <- all_positions[pos_idx]
    
    if (verbose && pos_idx %% 100 == 0) {
      cat("Processing position", pos_idx, "of", length(all_positions), "(", pos, ")\n")
    }
    
    # Subset data for this position ± 250kb (larger than max window)
    window_start <- max(1, pos - max_window_bp)
    window_end <- pos + max_window_bp
    
    # Find positions in window
    positions <- wide_data$refalt$POS
    window_positions <- positions >= window_start & positions <= window_end
    
    if (sum(window_positions) == 0) {
      if (verbose) cat("  No data in window for position", pos, "\n")
      next
    }
    
    # Create window data
    window_data <- list(
      refalt = wide_data$refalt[window_positions, ],
      freq_matrix = wide_data$freq_matrix[window_positions, , drop = FALSE],
      n_matrix = wide_data$n_matrix[window_positions, , drop = FALSE],
      sample_names = wide_data$sample_names
    )
    
    # Process all samples for this position
    for (sample_name in names_in_bam) {
      if (verbose && debug_level >= 2) {
        cat("  Processing sample:", sample_name, "at position:", pos, "\n")
      }
      
      # Call wide format haplotype estimation
      result <- estimate_haplotypes_wide(
        pos = pos,
        sample_name = sample_name,
        wide_data = window_data,
        founders = founders,
        h_cutoff = parameter,
        method = method,
        chr = chr,
        verbose = if (debug_level >= 2) 2 else 0
      )
      
      if (!is.null(result)) {
        result_counter <- result_counter + 1
        results_list[[result_counter]] <- list(
          CHROM = chr,
          pos = pos,
          sample = sample_name,
          Groups = result$Groups,
          Haps = result$Haps,
          Err = result$Err,
          Names = result$Names
        )
        
        if (verbose && debug_level >= 1) {
          cat("    Groups:", length(unique(result$Groups)), "founders distinguished\n")
          if (debug_level >= 2) {
            cat("    Haplotypes:", paste(names(result$Haps), "=", round(result$Haps, 3), collapse = ", "), "\n")
          }
        }
      } else {
        if (verbose && debug_level >= 1) {
          cat("    No result for sample:", sample_name, "at position:", pos, "\n")
        }
      }
    }
  }
  
  # Convert results to data frame
  if (length(results_list) == 0) {
    cat("No results generated!\n")
    return(data.frame())
  }
  
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      CHROM = x$CHROM,
      pos = x$pos,
      sample = x$sample,
      Groups = I(list(x$Groups)),
      Haps = I(list(x$Haps)),
      Err = I(list(x$Err)),
      Names = I(list(x$Names)),
      stringsAsFactors = FALSE
    )
  }))
  
  cat("Generated", nrow(results_df), "results\n")
  
  # Save results
  output_file <- file.path(output_dir, paste0("adaptive_window_h4_results_", chr, ".RDS"))
  saveRDS(results_df, output_file)
  cat("Results saved to:", output_file, "\n")
  
  return(results_df)
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
    cat("Limited to 100 positions × 1 sample for testing\n")
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
  
  adaptive_results <- run_adaptive_estimation_wide(chr, method, parameter, output_dir, param_file, debug, verbose, debug_level)
  
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
