#!/usr/bin/env Rscript

# Progressive Error Matrix Algorithm
# Implements the correct clustering threshold and founder distinguishability

library(limSolve)

# Function to create realistic founder states with proper distinguishability
create_realistic_founders <- function(n_obs_total, n_founders, h_cutoff = 4) {
  # Create founders with realistic distinguishability patterns
  A <- matrix(0, nrow = n_obs_total, ncol = n_founders)
  
  # F1: random vector
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
  
  return(A)
}

# Function to perform clustering with correct threshold
cluster_founders <- function(A_window, h_cutoff = 4) {
  founder_dist <- dist(t(A_window), method = "euclidean")
  hclust_result <- hclust(founder_dist, method = "complete")
  groups <- cutree(hclust_result, h = h_cutoff)  # Use h = h_cutoff, not h_cutoff^2
  return(groups)
}

# Function to build progressive error matrix
build_progressive_error_matrix <- function(A, window_sizes, h_cutoff = 4) {
  n_founders <- ncol(A)
  founder_names <- paste0("F", 1:n_founders)
  
  # Initialize error matrix
  error_matrix <- matrix(0, nrow = n_founders, ncol = n_founders)
  rownames(error_matrix) <- founder_names
  colnames(error_matrix) <- founder_names
  
  # Track which founders are locked in
  locked_founders <- rep(FALSE, n_founders)
  
  cat("=== PROGRESSIVE ERROR MATRIX BUILDING ===\n")
  
  for (window_size in window_sizes) {
    cat("\n--- Window size:", window_size, "SNPs ---\n")
    
    A_window <- A[1:window_size, ]
    groups <- cluster_founders(A_window, h_cutoff)
    
    # Show group composition
    unique_groups <- unique(groups)
    n_groups <- length(unique_groups)
    cat("Number of groups:", n_groups, "\n")
    
    for (i in 1:n_groups) {
      group_founders <- which(groups == unique_groups[i])
      group_name <- paste(founder_names[group_founders], collapse = "+")
      cat("  Group", i, ":", group_name, "(", length(group_founders), "founders)\n")
    }
    
    # Update error matrix for this window size
    for (i in 1:n_founders) {
      for (j in 1:n_founders) {
        if (i != j && !locked_founders[i] && !locked_founders[j]) {
          if (groups[i] == groups[j]) {
            # Founders are in same group - they are indistinguishable
            error_matrix[i, j] <- 0
          } else {
            # Founders are in different groups - they are distinguishable
            error_matrix[i, j] <- 1
          }
        }
      }
    }
    
    # Lock in founders that are now distinguishable from all others
    for (i in 1:n_founders) {
      if (!locked_founders[i]) {
        # Check if this founder is distinguishable from all others
        distinguishable <- TRUE
        for (j in 1:n_founders) {
          if (i != j && !locked_founders[j]) {
            if (groups[i] == groups[j]) {
              distinguishable <- FALSE
              break
            }
          }
        }
        if (distinguishable) {
          locked_founders[i] <- TRUE
          cat("  Locked in founder", founder_names[i], "\n")
        }
      }
    }
    
    # Show current error matrix
    cat("Current error matrix:\n")
    print(error_matrix)
  }
  
  return(error_matrix)
}

# Main function
main <- function() {
  # Set seed for reproducibility
  set.seed(123)
  
  # Simulation parameters
  n_obs_total <- 3000
  n_founders <- 8
  h_cutoff <- 4
  
  # Create realistic founders
  A <- create_realistic_founders(n_obs_total, n_founders, h_cutoff)
  
  # Define window sizes for progressive building
  window_sizes <- c(150, 300, 750, 1500, 3000)
  
  # Build progressive error matrix
  error_matrix <- build_progressive_error_matrix(A, window_sizes, h_cutoff)
  
  cat("\n=== FINAL ERROR MATRIX ===\n")
  print(error_matrix)
  
  # Check condition number
  condition_number <- kappa(error_matrix)
  cat("\nCondition number:", condition_number, "\n")
  
  if (condition_number < 1e12) {
    cat("✅ Error matrix is well-conditioned\n")
  } else {
    cat("❌ Error matrix is ill-conditioned\n")
  }
}

# Run if called directly
if (sys.nframe() == 0) {
  main()
}
