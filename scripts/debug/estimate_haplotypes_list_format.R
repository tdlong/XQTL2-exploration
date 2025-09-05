#!/usr/bin/env Rscript

# Modified estimate_haplotypes function that returns list format
# and shows adaptive algorithm progress

library(tidyverse)
library(limSolve)

#' Estimate Haplotype Frequencies with List Format Output
#' 
#' Modified version that returns list format and shows adaptive algorithm progress
#' 
#' @param pos Genomic position (integer)
#' @param sample_name Sample identifier (string)
#' @param df3 Processed data frame with POS, name, freq columns
#' @param founders Vector of founder names
#' @param h_cutoff Hierarchical clustering cutoff threshold
#' @param method Either "fixed" or "adaptive"
#' @param window_size_bp Window size in base pairs (required for fixed method)
#' @param chr Chromosome name (for output)
#' @param verbose Verbosity level (0=silent, 1=basic, 2=detailed, 3=full debug)
#' 
#' @return List with CHROM, pos, sample, Groups, Haps, Err, Names
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                                           method = c("fixed", "adaptive"),
                                           window_size_bp = NULL,
                                           chr = "chr2R",
                                           verbose = 0) {
  
  method <- match.arg(method)
  
  if (verbose >= 1) {
    cat("Processing pos:", pos, "sample:", sample_name, "method:", method, "\n")
  }
  
  # Initialize result variables
  groups <- rep(1, length(founders))
  names(groups) <- founders
  error_matrix <- matrix(NA, length(founders), length(founders))
  rownames(error_matrix) <- founders
  colnames(error_matrix) <- founders
  
  if (method == "fixed") {
    # Fixed window method
    if (verbose >= 1) {
      cat("Fixed window size:", window_size_bp, "bp\n")
    }
    
    # Get data for this position and sample
    sample_data <- df3 %>%
      dplyr::filter(POS >= pos - window_size_bp/2, POS <= pos + window_size_bp/2) %>%
      dplyr::filter(name == sample_name)
    
    if (nrow(sample_data) == 0) {
      if (verbose >= 1) {
        cat("No data found for position", pos, "sample", sample_name, "\n")
      }
      return(create_list_result(chr, pos, sample_name, method, window_size_bp, 0, NA, 
                               rep(NA, length(founders)), founders, groups, error_matrix))
    }
    
    # Run LSEI and clustering
    result <- run_lsei_and_clustering_list_format(sample_data, founders, h_cutoff, verbose)
    
    if (result$estimate_OK) {
      groups <- result$groups
      error_matrix <- result$error_matrix
    }
    
    return(create_list_result(chr, pos, sample_name, method, window_size_bp, 
                             nrow(sample_data), result$estimate_OK, result$haplotype_freqs, 
                             founders, groups, error_matrix))
    
  } else {
    # Adaptive window method
    if (verbose >= 1) {
      cat("=== ADAPTIVE WINDOW METHOD: h_cutoff =", h_cutoff, "===\n")
      cat("Position:", pos, "Sample:", sample_name, "\n")
    }
    
    # Define window sizes to test
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
    
    if (verbose >= 1) {
      cat("Window sizes to test:", paste(window_sizes/1000, "kb", collapse = " "), "\n")
    }
    
    final_result <- NULL
    final_window_size <- window_sizes[1]
    final_groups <- groups
    final_error_matrix <- error_matrix
    
    # Test each window size
    for (i in seq_along(window_sizes)) {
      window_size <- window_sizes[i]
      
      if (verbose >= 1) {
        cat("--- Window", i, ":", window_size/1000, "kb ---\n")
      }
      
      # Get data for this window size
      sample_data <- df3 %>%
        dplyr::filter(POS >= pos - window_size/2, POS <= pos + window_size/2) %>%
        dplyr::filter(name == sample_name)
      
      if (nrow(sample_data) < 10) {
        if (verbose >= 1) {
          cat("Insufficient SNPs (", nrow(sample_data), ") for window", window_size/1000, "kb\n")
        }
        next
      }
      
      # Run LSEI and clustering
      result <- run_lsei_and_clustering_list_format(sample_data, founders, h_cutoff, verbose)
      
      if (result$estimate_OK) {
        if (verbose >= 1) {
          cat("✓ Window", window_size/1000, "kb successful -", nrow(sample_data), "SNPs\n")
        }
        final_result <- result
        final_window_size <- window_size
        final_groups <- result$groups
        final_error_matrix <- result$error_matrix
        break
      } else {
        if (verbose >= 1) {
          cat("✗ Window", window_size/1000, "kb failed\n")
        }
      }
    }
    
    if (is.null(final_result)) {
      if (verbose >= 1) {
        cat("❌ All window sizes failed\n")
      }
      return(create_list_result(chr, pos, sample_name, method, final_window_size, 0, NA,
                               rep(NA, length(founders)), founders, final_groups, final_error_matrix))
    }
    
    if (verbose >= 1) {
      cat("✅ FINAL RESULT: estimate_OK:", final_result$estimate_OK, 
          "final_window_size:", final_window_size/1000, "kb n_snps:", nrow(sample_data), "\n")
    }
    
    return(create_list_result(chr, pos, sample_name, method, final_window_size,
                             nrow(sample_data), final_result$estimate_OK, final_result$haplotype_freqs,
                             founders, final_groups, final_error_matrix))
  }
}

#' Run LSEI estimation and clustering analysis with list format output
run_lsei_and_clustering_list_format <- function(sample_data, founders, h_cutoff, verbose) {
  
  # Create founder matrix and sample frequencies
  founder_matrix <- sample_data %>%
    dplyr::select(POS, all_of(founders)) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::select(-POS)
  
  sample_freqs <- sample_data %>%
    dplyr::select(POS, freq) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::pull(freq)
  
  if (nrow(founder_matrix) < 10) {
    return(list(estimate_OK = FALSE, haplotype_freqs = rep(NA, length(founders)), 
                groups = rep(1, length(founders)), error_matrix = matrix(NA, length(founders), length(founders))))
  }
  
  # Run LSEI with fulloutput=TRUE to capture error matrix
  result <- limSolve::lsei(A = as.matrix(founder_matrix), B = sample_freqs,
                          G = diag(length(founders)), H = matrix(rep(0.0003, length(founders))),
                          fulloutput = TRUE)
  
  if (result$IsError) {
    return(list(estimate_OK = FALSE, haplotype_freqs = rep(NA, length(founders)),
                groups = rep(1, length(founders)), error_matrix = matrix(NA, length(founders), length(founders))))
  }
  
  # Get haplotype frequencies
  haplotype_freqs <- result$X
  names(haplotype_freqs) <- founders
  
  # Normalize to sum to 1
  haplotype_freqs <- haplotype_freqs / sum(haplotype_freqs)
  
  # Run clustering
  distances <- dist(haplotype_freqs)
  hclust_result <- hclust(distances)
  groups <- cutree(hclust_result, h = h_cutoff)
  names(groups) <- founders
  
  # Get error matrix
  error_matrix <- result$cov
  if (!is.null(error_matrix)) {
    rownames(error_matrix) <- founders
    colnames(error_matrix) <- founders
  } else {
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders
    colnames(error_matrix) <- founders
  }
  
  # Check if founders are distinguishable
  n_groups <- length(unique(groups))
  estimate_OK <- (n_groups == length(founders))
  
  return(list(estimate_OK = estimate_OK, haplotype_freqs = haplotype_freqs,
              groups = groups, error_matrix = error_matrix))
}

#' Create result structure in list format
create_list_result <- function(chr, pos, sample_name, method, window_size, n_snps, estimate_OK, 
                              haplotype_freqs, founders, groups, error_matrix) {
  
  return(list(
    CHROM = chr,
    pos = pos,
    sample = sample_name,
    Groups = list(groups),           # Vector of length 8
    Haps = list(haplotype_freqs),    # Vector of length 8
    Err = list(error_matrix),        # Matrix 8x8
    Names = list(founders)           # Vector of length 8
  ))
}
