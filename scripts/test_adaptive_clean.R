#!/usr/bin/env Rscript

# Clean test script for adaptive window algorithm
# Shows: window position, window size, constraints, LSEI status

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

# Function to check if groups meaningfully changed
groups_changed <- function(current_groups, previous_groups) {
  if (is.null(previous_groups)) return(TRUE)
  
  # Check if number of groups changed
  if (length(unique(current_groups)) != length(unique(previous_groups))) return(TRUE)
  
  # Check if group composition changed
  current_group_list <- split(names(current_groups), current_groups)
  previous_group_list <- split(names(previous_groups), previous_groups)
  
  # Sort by group size and founder names for comparison
  current_sorted <- lapply(current_group_list, function(x) sort(x))
  previous_sorted <- lapply(previous_group_list, function(x) sort(x))
  
  # Compare sorted group compositions
  !identical(current_sorted, previous_sorted)
}

# Function to run adaptive window algorithm for one position and h_cutoff
run_adaptive_window_test <- function(test_pos, h_cutoff, sample_name) {
  cat("=== Testing position", test_pos, "with h_cutoff", h_cutoff, "for sample", sample_name, "===\n")
  
  # Initialize constraints for this h_cutoff
  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL
  
  # Track previous groups for comparison
  previous_groups <- NULL
  best_result <- NULL
  best_n_groups <- 0
  
  # Run adaptive window algorithm
  for (window_idx in seq_along(window_sizes)) {
    current_window_size <- window_sizes[window_idx]
    
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
      cat(sprintf("window position = %d | window size = %d kb | constraints = %s | lsei = %s\n", 
                  test_pos, current_window_size/1000, "insufficient SNPs", "skipped"))
      next
    }
    
    # Convert to matrix for clustering
    founder_matrix <- as.matrix(founder_matrix)
    
    # Cluster founders based on similarity using hierarchical clustering
    founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
    
    # Check if groups meaningfully changed
    groups_meaningfully_changed <- groups_changed(founder_clusters, previous_groups)
    
    # Show progress in condensed format
    n_groups <- length(unique(founder_clusters))
    constraint_info <- if (is.null(accumulated_constraints)) "none" else paste(nrow(accumulated_constraints), "constraints")
    
    if (groups_meaningfully_changed) {
      # Run LSEI with new grouping
      n_founders <- ncol(founder_matrix)
      E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
      F <- 1.0
      
      # Add accumulated constraints from previous windows
      if (!is.null(accumulated_constraints)) {
        E <- rbind(E, accumulated_constraints)
        F <- c(F, accumulated_constraint_values)
      }
      
      # Solve constrained least squares
      tryCatch({
        result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
        
        if (result$IsError == 0) {
          lsei_status <- "success"
          
          # Accumulate constraints for next window
          current_constraints <- NULL
          current_constraint_values <- NULL
          
          for (cluster_id in unique(founder_clusters)) {
            cluster_founders <- which(founder_clusters == cluster_id)
            if (length(cluster_founders) > 1) {
              # Group constraint
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              group_freq <- sum(result$X[cluster_founders])
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, group_freq)
            } else {
              # Individual constraint
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              founder_freq <- result$X[cluster_founders]
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
            }
          }
          
          # Update accumulated constraints
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
          }
          
          # Store best result
          best_result <- result
          best_n_groups <- n_groups
          
          # Check if all founders separated
          if (n_groups == length(founders)) {
            cat(sprintf("window position = %d | window size = %d kb | constraints = %s | lsei = %s | STOPPING (all founders separated)\n", 
                        test_pos, current_window_size/1000, constraint_info, lsei_status))
            break
          }
          
        } else {
          lsei_status <- "failed"
        }
      }, error = function(e) {
        lsei_status <- "error"
      })
      
    } else {
      # Groups unchanged, skip LSEI
      lsei_status <- "skipped (groups unchanged)"
    }
    
    # Show progress
    cat(sprintf("window position = %d | window size = %d kb | constraints = %s | lsei = %s\n", 
                test_pos, current_window_size/1000, constraint_info, lsei_status))
    
    # Update previous groups for next iteration
    previous_groups <- founder_clusters
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
