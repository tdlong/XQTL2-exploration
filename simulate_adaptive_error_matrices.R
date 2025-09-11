#!/usr/bin/env Rscript

# Comprehensive Simulation of Adaptive Method Variance-Covariance Matrices
# This script simulates the adaptive window method and demonstrates how error matrices
# are captured, averaged, and used in the smooth_h4 estimator

library(tidyverse)
library(limSolve)

cat("=== ADAPTIVE METHOD ERROR MATRIX SIMULATION ===\n\n")

# Set seed for reproducibility
set.seed(123)

# Simulation parameters
n_obs_total <- 3000
n_founders <- 8
founder_names <- paste0("F", 1:n_founders)
h_cutoff <- 4

# Create simulated founder haplotypes with different similarity levels
cat("1. Creating simulated founder haplotypes...\n")
A <- matrix(0, nrow = n_obs_total, ncol = n_founders)

# F1: random vector (baseline)
A[, 1] <- rbinom(n_obs_total, 1, 0.5)

# F2: 2 flips per 150 (distinguishable at 1500)
A[, 2] <- A[, 1]
block_size <- 150
flips_per_block <- 2
n_blocks <- n_obs_total %/% block_size
for (block in 1:n_blocks) {
  start_pos <- (block - 1) * block_size + 1
  end_pos <- min(block * block_size, n_obs_total)
  block_indices <- start_pos:end_pos
  flip_indices <- sample(block_indices, size = min(flips_per_block, length(block_indices)))
  A[flip_indices, 2] <- 1 - A[flip_indices, 2]
}

# F3: 10 flips per 150 (distinguishable at 300)
A[, 3] <- A[, 1]
flips_per_block <- 10
for (block in 1:n_blocks) {
  start_pos <- (block - 1) * block_size + 1
  end_pos <- min(block * block_size, n_obs_total)
  block_indices <- start_pos:end_pos
  flip_indices <- sample(block_indices, size = min(flips_per_block, length(block_indices)))
  A[flip_indices, 3] <- 1 - A[flip_indices, 3]
}

# F4: random vector
A[, 4] <- rbinom(n_obs_total, 1, 0.5)

# F5: 4 flips per 150 (distinguishable at 750)
A[, 5] <- A[, 4]
flips_per_block <- 4
for (block in 1:n_blocks) {
  start_pos <- (block - 1) * block_size + 1
  end_pos <- min(block * block_size, n_obs_total)
  block_indices <- start_pos:end_pos
  flip_indices <- sample(block_indices, size = min(flips_per_block, length(block_indices)))
  A[flip_indices, 5] <- 1 - A[flip_indices, 5]
}

# F6, F7, F8: random vectors
A[, 6] <- rbinom(n_obs_total, 1, 0.5)
A[, 7] <- rbinom(n_obs_total, 1, 0.5)
A[, 8] <- rbinom(n_obs_total, 1, 0.5)

cat("✓ Founder setup complete\n")
cat("  F1: random vector (baseline)\n")
cat("  F2: 2 flips per 150 (distinguishable at 1500)\n")
cat("  F3: 10 flips per 150 (distinguishable at 300)\n")
cat("  F4: random vector\n")
cat("  F5: 4 flips per 150 (distinguishable at 750)\n")
cat("  F6, F7, F8: random vectors\n\n")

# Function to simulate LSEI with error matrix capture
simulate_lsei_with_error <- function(founder_matrix, sample_freqs, h_cutoff, verbose = FALSE) {
  n_founders <- ncol(founder_matrix)
  
  # Calculate distances and perform clustering
  founder_dist <- dist(t(founder_matrix), method = "euclidean")
  hclust_result <- hclust(founder_dist, method = "complete")
  groups <- cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))
  
  if (verbose) {
    cat("  Groups:", paste(groups, collapse = ", "), "\n")
    cat("  Number of groups:", n_groups, "\n")
  }
  
  # Create constraint matrix E and vector F for grouped founders
  E <- matrix(0, nrow = n_groups, ncol = n_founders)
  F <- rep(1, n_groups)
  
  for (i in 1:n_groups) {
    group_founders <- which(groups == i)
    E[i, group_founders] <- 1
  }
  
  # Solve with LSEI and capture error matrix
  result <- limSolve::lsei(
    A = founder_matrix, 
    B = sample_freqs,
    E = E, 
    F = F,
    G = diag(n_founders), 
    H = matrix(rep(0.0003, n_founders)),
    fulloutput = TRUE  # This captures the error matrix
  )
  
  if (result$IsError == 0) {
    # Extract error matrix (covariance matrix)
    error_matrix <- result$cov
    rownames(error_matrix) <- founder_names
    colnames(error_matrix) <- founder_names
    
    return(list(
      success = TRUE,
      haplotype_freqs = result$X,
      groups = groups,
      error_matrix = error_matrix,
      n_groups = n_groups
    ))
  } else {
    return(list(
      success = FALSE,
      error_code = result$IsError,
      groups = groups,
      n_groups = n_groups
    ))
  }
}

# Test different window sizes and capture error matrices
cat("2. Testing adaptive method at different window sizes...\n")
window_sizes <- c(150, 300, 750, 1500, 3000)
results_list <- list()

for (window_size in window_sizes) {
  cat("\n--- Window size:", window_size, "SNPs ---\n")
  
  # Extract window data
  A_window <- A[1:window_size, ]
  
  # Create a simulated sample (mix of F1 and F4)
  sample_freqs <- 0.6 * A_window[, 1] + 0.4 * A_window[, 4] + rnorm(window_size, 0, 0.01)
  sample_freqs <- pmax(0, pmin(1, sample_freqs))  # Clamp to [0,1]
  
  # Run adaptive method
  result <- simulate_lsei_with_error(A_window, sample_freqs, h_cutoff, verbose = TRUE)
  
  if (result$success) {
    cat("✓ LSEI successful\n")
    cat("  Haplotype frequencies:", round(result$haplotype_freqs, 3), "\n")
    cat("  Sum of frequencies:", round(sum(result$haplotype_freqs), 3), "\n")
    
    # Store results
    results_list[[as.character(window_size)]] <- result
    
    # Show error matrix properties
    cat("  Error matrix properties:\n")
    cat("    Dimensions:", dim(result$error_matrix)[1], "x", dim(result$error_matrix)[2], "\n")
    cat("    Diagonal (variances):", round(diag(result$error_matrix), 4), "\n")
    cat("    Trace (total variance):", round(sum(diag(result$error_matrix)), 4), "\n")
    cat("    Determinant:", round(det(result$error_matrix), 6), "\n")
    
    # Show off-diagonal elements (covariances)
    off_diag <- result$error_matrix
    diag(off_diag) <- NA
    cat("    Max off-diagonal (covariance):", round(max(off_diag, na.rm = TRUE), 4), "\n")
    cat("    Min off-diagonal (covariance):", round(min(off_diag, na.rm = TRUE), 4), "\n")
    
  } else {
    cat("✗ LSEI failed with error code:", result$error_code, "\n")
  }
}

# Demonstrate error matrix averaging (like in smooth_h4)
cat("\n3. Demonstrating error matrix averaging (smooth_h4 style)...\n")

# Get successful results
successful_results <- results_list[sapply(results_list, function(x) x$success)]

if (length(successful_results) > 0) {
  cat("✓ Found", length(successful_results), "successful results for averaging\n")
  
  # Extract error matrices
  error_matrices <- lapply(successful_results, function(x) x$error_matrix)
  
  # Calculate average error matrix
  avg_error_matrix <- Reduce(`+`, error_matrices) / length(error_matrices)
  
  cat("Average error matrix properties:\n")
  cat("  Dimensions:", dim(avg_error_matrix)[1], "x", dim(avg_error_matrix)[2], "\n")
  cat("  Diagonal (avg variances):", round(diag(avg_error_matrix), 4), "\n")
  cat("  Trace (total avg variance):", round(sum(diag(avg_error_matrix)), 4), "\n")
  cat("  Determinant:", round(det(avg_error_matrix), 6), "\n")
  
  # Show the full average error matrix
  cat("\nAverage error matrix:\n")
  print(round(avg_error_matrix, 4))
  
  # Calculate variance of variances across window sizes
  variances_by_window <- sapply(error_matrices, diag)
  if (ncol(variances_by_window) > 1) {
    cat("\nVariance of variances across window sizes:\n")
    for (i in 1:n_founders) {
      cat("  F", i, " variance range: [", 
          round(min(variances_by_window[i, ]), 4), ", ", 
          round(max(variances_by_window[i, ]), 4), "]\n")
    }
  }
  
} else {
  cat("✗ No successful results found for averaging\n")
}

# Demonstrate the 4-list structure for smooth_h4
cat("\n4. Demonstrating 4-list structure for smooth_h4...\n")

if (length(successful_results) > 0) {
  # Extract the last successful result (largest window with all founders distinguishable)
  final_result <- successful_results[[length(successful_results)]]
  
  # Create the 4 lists as used in smooth_h4
  groups_list <- final_result$groups
  haps_list <- final_result$haplotype_freqs
  err_list <- final_result$error_matrix
  names_list <- founder_names
  
  cat("✓ 4-list structure created:\n")
  cat("  1. Groups:", paste(groups_list, collapse = ", "), "\n")
  cat("  2. Haps:", round(haps_list, 4), "\n")
  cat("  3. Err dimensions:", dim(err_list)[1], "x", dim(err_list)[2], "\n")
  cat("  4. Names:", paste(names_list, collapse = ", "), "\n")
  
  # Show how this would be used in a tibble
  cat("\nTibble structure example:\n")
  example_tibble <- tibble(
    CHROM = "chr2R",
    pos = 5400000,
    sample = "Rep01_W_F",
    Groups = list(groups_list),
    Haps = list(haps_list),
    Err = list(err_list),
    Names = list(names_list)
  )
  
  print(example_tibble)
}

cat("\n=== SIMULATION COMPLETE ===\n")
cat("This simulation demonstrates:\n")
cat("1. How the adaptive method groups founders at different window sizes\n")
cat("2. How error matrices are captured from LSEI with fulloutput=TRUE\n")
cat("3. How error matrices are averaged across positions (smooth_h4 style)\n")
cat("4. The 4-list structure used in list-format haplotype estimation\n")

