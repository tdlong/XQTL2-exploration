#!/usr/bin/env Rscript

# Progressive Error Matrix Algorithm - Detailed Implementation
# Shows progressive fitting, error matrix building, and parameter estimates

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
  groups <- cutree(hclust_result, h = h_cutoff)
  return(groups)
}

# Function to fit LSEI and get parameter estimates
fit_lsei <- function(A_window, groups, founder_names) {
  n_founders <- ncol(A_window)
  n_obs <- nrow(A_window)
  
  # Create constraint matrix
  n_groups <- length(unique(groups))
  n_constraints <- n_groups
  
  # Initialize constraint matrix
  G <- matrix(0, nrow = n_constraints, ncol = n_founders)
  h <- rep(1, n_constraints)
  
  # Add constraints for each group
  for (i in 1:n_groups) {
    group_founders <- which(groups == i)
    G[i, group_founders] <- 1
  }
  
  # Create observation matrix (A_window)
  # Each row is an observation, each column is a founder
  # We want to estimate the proportion of each founder in each observation
  
  # For simplicity, let's create a target vector (this would be the actual data in practice)
  # Here we'll use the mean of each founder as the target
  target <- colMeans(A_window)
  
  # Fit LSEI
  tryCatch({
    result <- lsei(A = A_window, B = target, G = G, H = h, type = 2)
    
    # Extract parameter estimates
    params <- result$X
    names(params) <- founder_names
    
    # Extract error matrix
    error_matrix <- result$covar
    rownames(error_matrix) <- founder_names
    colnames(error_matrix) <- founder_names
    
    return(list(
      params = params,
      error_matrix = error_matrix,
      success = TRUE,
      message = "LSEI fit successful"
    ))
  }, error = function(e) {
    return(list(
      params = rep(NA, n_founders),
      error_matrix = matrix(NA, nrow = n_founders, ncol = n_founders),
      success = FALSE,
      message = paste("LSEI fit failed:", e$message)
    ))
  })
}

# Function to build progressive error matrix with detailed output
build_progressive_error_matrix_detailed <- function(A, window_sizes, h_cutoff = 4) {
  n_founders <- ncol(A)
  founder_names <- paste0("F", 1:n_founders)
  
  # Initialize error matrix
  error_matrix <- matrix(0, nrow = n_founders, ncol = n_founders)
  rownames(error_matrix) <- founder_names
  colnames(error_matrix) <- founder_names
  
  # Track which founders are locked in
  locked_founders <- rep(FALSE, n_founders)
  
  # Store results for each window size
  results <- list()
  
  cat("=== PROGRESSIVE ERROR MATRIX BUILDING - DETAILED ===\n")
  
  for (window_idx in 1:length(window_sizes)) {
    window_size <- window_sizes[window_idx]
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("WINDOW SIZE:", window_size, "SNPs\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    
    A_window <- A[1:window_size, ]
    groups <- cluster_founders(A_window, h_cutoff)
    
    # Show group composition
    unique_groups <- unique(groups)
    n_groups <- length(unique_groups)
    cat("\n--- CLUSTERING RESULTS ---\n")
    cat("Number of groups:", n_groups, "\n")
    
    for (i in 1:n_groups) {
      group_founders <- which(groups == i)
      group_name <- paste(founder_names[group_founders], collapse = "+")
      cat("  Group", i, ":", group_name, "(", length(group_founders), "founders)\n")
    }
    
    # Show distances
    founder_dist <- dist(t(A_window), method = "euclidean")
    dist_matrix <- as.matrix(founder_dist)
    cat("\n--- DISTANCE MATRIX ---\n")
    print(round(dist_matrix, 2))
    
    # Fit LSEI
    cat("\n--- LSEI FITTING ---\n")
    lsei_result <- fit_lsei(A_window, groups, founder_names)
    
    if (lsei_result$success) {
      cat("✅ LSEI fit successful\n")
      cat("Parameter estimates:\n")
      print(round(lsei_result$params, 4))
      
      cat("\nError matrix (covariance):\n")
      print(round(lsei_result$error_matrix, 4))
      
      # Check condition number
      condition_number <- kappa(lsei_result$error_matrix)
      cat("Condition number:", condition_number, "\n")
      
      if (condition_number < 1e12) {
        cat("✅ Error matrix is well-conditioned\n")
      } else {
        cat("❌ Error matrix is ill-conditioned\n")
      }
    } else {
      cat("❌ LSEI fit failed:", lsei_result$message, "\n")
    }
    
    # Update error matrix for this window size
    cat("\n--- ERROR MATRIX UPDATE ---\n")
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
    cat("\n--- FOUNDER LOCKING ---\n")
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
          cat("  🔒 Locked in founder", founder_names[i], "\n")
        }
      }
    }
    
    # Show current error matrix
    cat("\n--- CURRENT ERROR MATRIX ---\n")
    print(error_matrix)
    
    # Store results
    results[[window_idx]] <- list(
      window_size = window_size,
      groups = groups,
      n_groups = n_groups,
      lsei_result = lsei_result,
      error_matrix = error_matrix,
      locked_founders = locked_founders
    )
  }
  
  return(list(
    final_error_matrix = error_matrix,
    results = results
  ))
}

# Main function
main <- function() {
  # Set seed for reproducibility
  set.seed(123)
  
  # Simulation parameters
  n_obs_total <- 3000
  n_founders <- 8
  h_cutoff <- 4
  
  cat("=== PROGRESSIVE ERROR MATRIX ALGORITHM - DETAILED ===\n")
  cat("Simulation parameters:\n")
  cat("  Total observations:", n_obs_total, "\n")
  cat("  Number of founders:", n_founders, "\n")
  cat("  H cutoff:", h_cutoff, "\n")
  cat("  Window sizes: 150, 300, 750, 1500, 3000\n\n")
  
  # Create realistic founders
  cat("Creating realistic founders...\n")
  A <- create_realistic_founders(n_obs_total, n_founders, h_cutoff)
  
  # Define window sizes for progressive building
  window_sizes <- c(150, 300, 750, 1500, 3000)
  
  # Build progressive error matrix
  cat("Building progressive error matrix...\n")
  result <- build_progressive_error_matrix_detailed(A, window_sizes, h_cutoff)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("FINAL RESULTS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  cat("\n--- FINAL ERROR MATRIX ---\n")
  print(result$final_error_matrix)
  
  # Check final condition number
  condition_number <- kappa(result$final_error_matrix)
  cat("\nFinal condition number:", condition_number, "\n")
  
  if (condition_number < 1e12) {
    cat("✅ Final error matrix is well-conditioned\n")
  } else {
    cat("❌ Final error matrix is ill-conditioned\n")
  }
  
  # Show progression summary
  cat("\n--- PROGRESSION SUMMARY ---\n")
  for (i in 1:length(result$results)) {
    res <- result$results[[i]]
    cat("Window", res$window_size, "SNPs: ", res$n_groups, "groups,", 
        sum(res$locked_founders), "founders locked\n")
  }
}

# Run if called directly
if (sys.nframe() == 0) {
  main()
}
