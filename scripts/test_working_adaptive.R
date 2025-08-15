#!/usr/bin/env Rscript

# Test script with the working adaptive window algorithm from August 12th
# This tests the hierarchical clustering approach on a couple of positions

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Load parameter file
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")

# Load data
cat("Loading REFALT data...\n")
df <- read.table("process/JUICE/RefAlt.chr2R.txt", header = TRUE)

# Transform to frequencies
cat("Converting counts to frequencies...\n")
df2 <- df %>%
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

# Filter for high-quality SNPs
cat("Filtering for high-quality SNPs...\n")
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Subset dataset
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Non-founder samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Test positions in the middle of the chromosome
test_positions <- c(15000000, 20000000)  # 2 test positions
test_h_cutoffs <- c(4, 10)  # Test two different cutoffs

cat("Testing positions:", paste(test_positions, collapse = ", "), "\n")
cat("Testing h_cutoff values:", paste(test_h_cutoffs, collapse = ", "), "\n\n")

# Define window sizes to try (progressive expansion)
base_window_size <- 10000  # Start with 10kb
window_sizes <- c(base_window_size, 
                  base_window_size * 2.5, 
                  base_window_size * 5, 
                  base_window_size * 10, 
                  base_window_size * 20,
                  base_window_size * 50)  # Max 500kb

# Function to run adaptive window algorithm for one position and h_cutoff
run_adaptive_window_test <- function(test_pos, h_cutoff, sample_name) {
  cat("=== Testing position", test_pos, "with h_cutoff", h_cutoff, "for sample", sample_name, "===\n")
  
  # Initialize constraints for this h_cutoff
  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL
  
  # Track clustering progress
  previous_n_groups <- 0
  best_result <- NULL
  best_n_groups <- 0
  
  # Run adaptive window algorithm
  for (window_idx in seq_along(window_sizes)) {
    current_window_size <- window_sizes[window_idx]
    
    cat("  Window size:", current_window_size/1000, "kb\n")
    
    # Create window around test position
    window_start <- max(0, test_pos - current_window_size/2)
    window_end <- test_pos + current_window_size/2
    
    # Get SNPs in this expanding window
    window_snps <- df3 %>%
      filter(CHROM == "chr2R" &
             POS > window_start &
             POS < window_end &
             (name %in% founders | name == sample_name)) %>%
      select(-c(CHROM, N)) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    # Get founder matrix and sample frequencies
    founder_matrix <- window_snps %>% select(all_of(founders))
    sample_freqs <- window_snps[[sample_name]]
    
    # Filter for non-NA values
    valid_positions <- !is.na(sample_freqs)
    sample_freqs <- sample_freqs[valid_positions]
    founder_matrix <- founder_matrix[valid_positions, ]
    
    if (nrow(founder_matrix) < 10) {
      cat("    ⚠️  Insufficient SNPs in window\n")
      next
    }
    
    # Convert to matrix for clustering
    founder_matrix <- as.matrix(founder_matrix)
    
    # Cluster founders based on similarity using hierarchical clustering
    founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
    
    # Count groups
    n_groups <- length(unique(founder_clusters))
    
    cat("    Founder groups:", n_groups, "\n")
    
    # Check if clustering improved (more groups = better separation)
    if (window_idx > 1 && n_groups <= previous_n_groups) {
      cat("    Clustering not improved, skipping to next window\n")
      next
    }
    
    # Update previous_n_groups for next iteration
    previous_n_groups <- n_groups
    
    # Show group composition
    for (group_id in unique(founder_clusters)) {
      group_founders <- names(founder_clusters[founder_clusters == group_id])
      cat("      Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
    }
    
    # Build constraint matrix with accumulated constraints from smaller windows
    n_founders <- ncol(founder_matrix)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add accumulated constraints from previous (smaller) windows
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      cat("    Added", nrow(accumulated_constraints), "accumulated constraints\n")
    }
    
    # Get unique clusters
    unique_clusters <- unique(founder_clusters)
    
    # Solve constrained least squares
    tryCatch({
      result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                              G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
      
      if (result$IsError == 0) {
        cat("    ✓ LSEI successful\n")
        
        # Accumulate constraints for next (larger) window
        current_constraints <- NULL
        current_constraint_values <- NULL
        
        for (cluster_id in unique_clusters) {
          cluster_founders <- which(founder_clusters == cluster_id)
          if (length(cluster_founders) > 1) {
            # Create constraint row for this group
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            # Calculate the actual group frequency from lsei result
            group_freq <- sum(result$X[cluster_founders])
            
            cat("      Group", cluster_id, "frequency:", round(group_freq, 4), "\n")
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, group_freq)
          } else {
            # Single founder: lock their exact frequency
            founder_freq <- result$X[cluster_founders]
            
            cat("      Single founder frequency:", round(founder_freq, 4), "\n")
            
            # Create constraint: this founder = their exact frequency
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, founder_freq)
          }
        }
        
        # Update accumulated constraints for next window
        if (!is.null(current_constraints)) {
          accumulated_constraints <- current_constraints
          accumulated_constraint_values <- current_constraint_values
          cat("    ✓ Accumulated", nrow(current_constraints), "constraints for next window\n")
        } else {
          # If no constraints (e.g., all founders in 1 group), reset to NULL
          accumulated_constraints <- NULL
          accumulated_constraint_values <- NULL
          cat("    ✓ No meaningful constraints to accumulate\n")
        }
        
        # Store the best result for this position/h_cutoff
        best_result <- result
        best_n_groups <- n_groups
        
        # Check if all founders are separated
        if (n_groups == length(founders)) {
          cat("    ✓ All founders separated! Stopping for this h_cutoff.\n")
          break
        }
        
      } else {
        cat("    ❌ LSEI failed\n")
      }
    }, error = function(e) {
      cat("    ❌ LSEI error:", e$message, "\n")
    })
  }
  
  # Return the best result found
  if (!is.null(best_result)) {
    cat("  Final result: Found", best_n_groups, "founder groups\n")
    cat("  Founder frequencies:", round(best_result$X, 4), "\n")
    return(list(
      h_cutoff = h_cutoff,
      n_groups = best_n_groups,
      founder_frequencies = best_result$X,
      converged = TRUE
    ))
  } else {
    cat("  Final result: No successful estimation\n")
    return(list(
      h_cutoff = h_cutoff,
      n_groups = 0,
      founder_frequencies = rep(NA, length(founders)),
      converged = FALSE
    ))
  }
}

# Test each combination
results <- list()

for (test_pos in test_positions) {
  for (h_cutoff in test_h_cutoffs) {
    for (sample_name in head(non_founder_samples, 1)) {  # Test just first sample
      result_key <- paste0("pos_", test_pos, "_h", h_cutoff, "_", sample_name)
      results[[result_key]] <- run_adaptive_window_test(test_pos, h_cutoff, sample_name)
      cat("\n")
    }
  }
}

# Compare results
cat("=== COMPARISON OF RESULTS ===\n")

# Compare h4 vs h10 for each position
for (test_pos in test_positions) {
  cat("\nPosition", test_pos, ":\n")
  
  h4_key <- paste0("pos_", test_pos, "_h4_", head(non_founder_samples, 1))
  h10_key <- paste0("pos_", test_pos, "_h10_", head(non_founder_samples, 1))
  
  if (h4_key %in% names(results) && h10_key %in% names(results)) {
    h4_result <- results[[h4_key]]
    h10_result <- results[[h10_key]]
    
    cat("  h4:  n_groups =", h4_result$n_groups, "frequencies =", round(h4_result$founder_frequencies, 4), "\n")
    cat("  h10: n_groups =", h10_result$n_groups, "frequencies =", round(h10_result$founder_frequencies, 4), "\n")
    
    # Check if frequencies are identical
    if (h4_result$converged && h10_result$converged) {
      freq_diff <- max(abs(h4_result$founder_frequencies - h10_result$founder_frequencies), na.rm = TRUE)
      cat("  Max frequency difference:", round(freq_diff, 6), "\n")
      
      if (freq_diff < 1e-6) {
        cat("  ⚠️  IDENTICAL RESULTS - Algorithm may not be working properly\n")
      } else {
        cat("  ✓ DIFFERENT RESULTS - Algorithm is working!\n")
      }
    } else {
      cat("  ⚠️  One or both failed to converge\n")
    }
  }
}

cat("\n=== TEST COMPLETE ===\n")
