#!/usr/bin/env Rscript

# Centromere-Specific Haplotype Estimation Function
# This function implements a single-window approach for centromere regions
# It uses the full range of centromere SNPs as one window (no sliding/adaptive)

library(tidyverse)
library(limSolve)

#' Estimate Haplotype Frequencies for Centromere Regions
#' 
#' Single-window haplotype estimation for centromere regions using all available SNPs
#' 
#' @param pos Genomic position (integer) - center position for the window
#' @param sample_name Sample identifier (string)
#' @param df3 Processed data frame with POS, name, freq columns (already filtered to centromere SNPs)
#' @param founders Vector of founder names
#' @param h_cutoff Hierarchical clustering cutoff threshold
#' @param window_size_bp Window size in base pairs (uses all centromere SNPs)
#' @param chr Chromosome name (for output)
#' @param verbose Verbosity level (0=silent, 1=basic, 2=detailed, 3=full debug)
#' 
#' @return List with:
#'   - Groups: vector of group assignments from clustering (length = number of founders)
#'   - Haps: vector of founder haplotype frequencies (length = number of founders, sums to 1.0)
#'   - Err: error covariance matrix from lsei (n_founders x n_founders)
#'   - Names: vector of founder names (length = number of founders)
estimate_haplotypes_single_window <- function(pos, sample_name, df3, founders, h_cutoff,
                                          window_size_bp = NULL,
                                          chr = "chr2R",
                                          verbose = 0) {
  
  if (verbose >= 1) {
    cat(sprintf("CENTROMERE: Processing pos: %s, sample: %s\n", 
                format(pos, big.mark=","), sample_name))
  }
  
  # For centromere, we use ALL SNPs in df3 (already filtered to centromere region)
  # No window size needed - we use everything
  
  if (verbose >= 2) {
    cat(sprintf("=== CENTROMERE SINGLE WINDOW: h_cutoff = %g ===\n", h_cutoff))
    cat(sprintf("Position: %s\n", format(pos, big.mark=",")))
    cat(sprintf("Sample: %s\n", sample_name))
    cat(sprintf("Using ALL centromere SNPs: %d positions\n", length(unique(df3$POS))))
  }
  
  # Get all positions in the centromere region
  centromere_positions <- unique(df3$POS)
  
  if (length(centromere_positions) == 0) {
    if (verbose >= 1) {
      cat("ERROR: No centromere positions found\n")
    }
    return(NULL)
  }
  
  # Filter data to centromere positions and current sample
  sample_data <- df3 %>%
    filter(POS %in% centromere_positions, name == sample_name) %>%
    select(POS, freq) %>%
    arrange(POS)
  
  if (nrow(sample_data) == 0) {
    if (verbose >= 1) {
      cat("ERROR: No data for sample", sample_name, "in centromere region\n")
    }
    return(NULL)
  }
  
  # Get founder data for all centromere positions
  founder_data <- df3 %>%
    filter(POS %in% centromere_positions, name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq) %>%
    arrange(POS)
  
  if (nrow(founder_data) == 0) {
    if (verbose >= 1) {
      cat("ERROR: No founder data in centromere region\n")
    }
    return(NULL)
  }
  
  # Check if we have enough data
  if (nrow(founder_data) < 3) {
    if (verbose >= 1) {
      cat("ERROR: Not enough centromere positions (", nrow(founder_data), ")\n")
    }
    return(NULL)
  }
  
  if (verbose >= 2) {
    cat("Centromere data summary:\n")
    cat("  Positions:", nrow(founder_data), "\n")
    cat("  Founders:", length(founders), "\n")
    cat("  Sample data points:", nrow(sample_data), "\n")
  }
  
  # Prepare data for haplotype estimation
  # X matrix: founder frequencies at each position
  X <- as.matrix(founder_data[, founders])
  
  # y vector: sample frequencies at each position
  y <- sample_data$freq[match(founder_data$POS, sample_data$POS)]
  
  # Remove any positions with missing data
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  if (length(y) < 3) {
    if (verbose >= 1) {
      cat("ERROR: Not enough complete positions after filtering (", length(y), ")\n")
    }
    return(NULL)
  }
  
  if (verbose >= 2) {
    cat("After filtering: ", length(y), "complete positions\n")
  }
  
  # Hierarchical clustering on founders (not positions)
  # Calculate distance matrix between founders based on their frequencies across positions
  founder_dist <- dist(t(X), method = "euclidean")
  
  # Print distance matrix for debugging
  cat("Founder distance matrix:\n")
  dist_matrix <- as.matrix(founder_dist)
  rownames(dist_matrix) <- founders
  colnames(dist_matrix) <- founders
  print(round(dist_matrix, 4))
  cat("\n")
  
  # Check correlation between B1 and AB8 specifically
  if ("B1" %in% founders && "AB8" %in% founders) {
    b1_idx <- which(founders == "B1")
    ab8_idx <- which(founders == "AB8")
    correlation <- cor(X[, b1_idx], X[, ab8_idx])
    cat("B1 vs AB8 correlation:", round(correlation, 4), "\n")
    cat("B1 + AB8 sum range:", round(range(X[, b1_idx] + X[, ab8_idx]), 4), "\n")
    cat("B1 range:", round(range(X[, b1_idx]), 4), "\n")
    cat("AB8 range:", round(range(X[, ab8_idx]), 4), "\n\n")
  }
  
  hc <- hclust(founder_dist, method = "complete")
  groups <- cutree(hc, h = h_cutoff)
  
  cat("Clustering result: ", length(unique(groups)), "groups\n")
  cat("Group assignments:", paste(groups, collapse = " "), "\n")
  cat("h_cutoff used:", h_cutoff, "\n")
  
  # Show which founders are in each group
  for (group_id in unique(groups)) {
    group_founders <- founders[groups == group_id]
    cat("Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
  }
  cat("\n")
  
  # Solve for haplotype frequencies using lsei
  # Constraint 1: frequencies must sum to 1
  E <- matrix(1, nrow = 1, ncol = length(founders))
  f <- 1
  
  # Constraint 2: frequencies must be non-negative (>= 0.0003)
  G <- diag(length(founders))
  h <- rep(0.0003, length(founders))
  
  # Solve the constrained least squares problem
  result <- lsei(A = X, B = y, E = E, F = f, G = G, H = h, fulloutput = TRUE)
  
  if (result$IsError) {
    if (verbose >= 1) {
      cat("ERROR: lsei failed to solve\n")
    }
    return(NULL)
  }
  
  # Extract founder frequencies
  founder_frequencies <- result$X
  names(founder_frequencies) <- founders
  
  # Normalize to ensure they sum to 1
  founder_frequencies <- founder_frequencies / sum(founder_frequencies)
  
  if (verbose >= 2) {
    cat("Haplotype frequencies:\n")
    for (i in seq_along(founders)) {
      cat(sprintf("  %s: %.4f\n", founders[i], founder_frequencies[i]))
    }
    cat("Sum:", sum(founder_frequencies), "\n")
  }
  
  # Return results in the expected format
  return(list(
    Groups = groups,
    Haps = founder_frequencies,
    Err = result$cov,
    Names = founders
  ))
}
