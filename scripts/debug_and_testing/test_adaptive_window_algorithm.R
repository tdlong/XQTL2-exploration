#!/usr/bin/env Rscript

# Complete Test Script - Adaptive Window Algorithm
# Tests the FULL adaptive window algorithm with progressive expansion and constraint accumulation

library(tidyverse)
library(limSolve)

cat("=== COMPLETE ADAPTIVE WINDOW ALGORITHM TEST ===\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and transform data (same as production)
cat("\n2. Loading and transforming data...\n")
filein <- "process/JUICE/RefAlt.chr2R.txt"
df <- read.table(filein, header = TRUE)

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
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("✓ Data ready:", nrow(df3), "rows\n")

# 3. Test the complete adaptive window algorithm
cat("\n3. Testing complete adaptive window algorithm...\n")

# Test parameters
test_pos <- 15000000
test_h_cutoff <- 4
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test h_cutoff:", test_h_cutoff, "\n")
cat("✓ Test sample:", test_sample, "\n\n")

# Define progressive window sizes (same as production)
base_window_size <- 10000
window_sizes <- c(base_window_size, 
                  base_window_size * 2.5, 
                  base_window_size * 5, 
                  base_window_size * 10, 
                  base_window_size * 20,
                  base_window_size * 50)

cat("Window sizes to test:", paste(window_sizes/1000, "kb"), "\n\n")

# Function to check if groups meaningfully changed (same as production)
groups_changed <- function(current_groups, previous_groups) {
  if (is.null(previous_groups)) return(TRUE)
  if (length(unique(current_groups)) != length(unique(previous_groups))) return(TRUE)
  current_group_list <- split(names(current_groups), current_groups)
  previous_group_list <- split(names(previous_groups), previous_groups)
  current_sorted <- lapply(current_group_list, function(x) sort(x))
  previous_sorted <- lapply(previous_group_list, function(x) sort(x))
  !identical(current_sorted, previous_sorted)
}

# Initialize for this h_cutoff
accumulated_constraints <- NULL
accumulated_constraint_values <- NULL
previous_groups <- NULL
best_result <- NULL
best_n_groups <- 0

cat("=== STARTING ADAPTIVE WINDOW ALGORITHM ===\n")

# Run the complete adaptive window algorithm
for (window_idx in seq_along(window_sizes)) {
  current_window_size <- window_sizes[window_idx]
  
  cat(sprintf("\n--- Window %d: %d kb ---\n", window_idx, current_window_size/1000))
  
  # Create window around test position
  window_start <- max(0, test_pos - current_window_size/2)
  window_end <- test_pos + current_window_size/2
  
  cat("Window range:", window_start, "to", window_end, "bp\n")
  
  # Get SNPs in this expanding window (same as production script)
  window_snps <- df3 %>%
    filter(CHROM == "chr2R" &
           POS > window_start &
           POS < window_end &
           (name %in% founders | name == test_sample))
  
  cat("SNPs in window:", nrow(window_snps), "\n")
  
  # Get founder data for this window (same as production script)
  founder_data <- window_snps %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  # Check if we have enough founder data
  if (ncol(founder_data) < length(founders) + 1) {
    cat("✗ Insufficient founder data, skipping\n")
    next
  }
  
  # Get founder matrix and sample frequencies
  founder_matrix <- founder_data %>%
    select(-POS) %>%
    as.matrix()
  
  # Get sample frequencies for the same positions
  sample_freqs <- window_snps %>%
    filter(name == test_sample) %>%
    select(POS, freq) %>%
    right_join(founder_data %>% select(POS), by = "POS") %>%
    pull(freq)
  
  # Filter for non-NA values
  valid_positions <- !is.na(sample_freqs)
  sample_freqs <- sample_freqs[valid_positions]
  founder_matrix <- founder_matrix[valid_positions, ]
  
  if (nrow(founder_matrix) < 10) {
    cat("✗ Insufficient data for estimation, skipping\n")
    next
  }
  
  cat("Valid positions for estimation:", nrow(founder_matrix), "\n")
  
  # Convert to matrix for clustering
  founder_matrix <- as.matrix(founder_matrix)
  
  # Cluster founders based on similarity using hierarchical clustering
  founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = test_h_cutoff)
  n_groups <- length(unique(founder_clusters))
  
  cat("Founder groups:", n_groups, "\n")
  cat("Group assignments:", paste(founder_clusters, collapse=", "), "\n")
  
  # Check if groups meaningfully changed
  groups_meaningfully_changed <- groups_changed(founder_clusters, previous_groups)
  cat("Groups meaningfully changed:", groups_meaningfully_changed, "\n")
  
  if (groups_meaningfully_changed) {
    cat("Running LSEI with new grouping...\n")
    
    # Run LSEI with new grouping
    n_founders <- ncol(founder_matrix)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add accumulated constraints from previous windows
    if (!is.null(accumulated_constraints)) {
      cat("Adding", nrow(accumulated_constraints), "accumulated constraints\n")
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
    }
    
    cat("Constraint matrix E:", nrow(E), "x", ncol(E), "\n")
    
    # Solve constrained least squares
    tryCatch({
      result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                              G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
      
      if (result$IsError == 0) {
        cat("✓ LSEI successful!\n")
        cat("Estimated frequencies:", paste(round(result$X, 4), collapse=", "), "\n")
        cat("Sum of frequencies:", sum(result$X), "\n")
        
        # Accumulate constraints for next (larger) window
        current_constraints <- NULL
        current_constraint_values <- NULL
        
        for (cluster_id in unique(founder_clusters)) {
          cluster_founders <- which(founder_clusters == cluster_id)
          if (length(cluster_founders) > 1) {
            # Create constraint row for this group
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            # Calculate the actual group frequency from lsei result
            group_freq <- sum(result$X[cluster_founders])
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, group_freq)
            
            cat("Group", cluster_id, "constraint:", paste(founders[cluster_founders], collapse="+"), "=", round(group_freq, 4), "\n")
          } else {
            # Single founder: lock their exact frequency
            founder_freq <- result$X[cluster_founders]
            
            # Create constraint: this founder = their exact frequency
            constraint_row <- rep(0, n_founders)
            constraint_row[cluster_founders] <- 1
            
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, founder_freq)
            
            cat("Founder", founders[cluster_founders], "locked at:", round(founder_freq, 4), "\n")
          }
        }
        
        # Update accumulated constraints for next window
        if (!is.null(current_constraints)) {
          accumulated_constraints <- current_constraints
          accumulated_constraint_values <- current_constraint_values
          cat("✓ Accumulated", nrow(accumulated_constraints), "constraints for next window\n")
        } else {
          # If no constraints (e.g., all founders in 1 group), reset to NULL
          accumulated_constraints <- NULL
          accumulated_constraint_values <- NULL
          cat("Reset constraints (all founders in one group)\n")
        }
        
        # Store the best result for this position/h_cutoff
        best_result <- result
        best_n_groups <- length(unique(founder_clusters))
        
        # Check if all founders are separated
        if (best_n_groups == length(founders)) {
          cat("✓ All founders separated! Stopping for this h_cutoff\n")
          break
        }
      } else {
        cat("✗ LSEI failed with error code:", result$IsError, "\n")
      }
    }, error = function(e) {
      cat("✗ LSEI error:", e$message, "\n")
    })
  } else {
    cat("Groups unchanged, skipping LSEI\n")
  }
  
  previous_groups <- founder_clusters
}

# Final results
cat("\n=== ADAPTIVE WINDOW ALGORITHM COMPLETE ===\n")
if (!is.null(best_result)) {
  cat("✓ Algorithm completed successfully\n")
  cat("Final number of groups:", best_n_groups, "\n")
  cat("Final estimated frequencies:", paste(round(best_result$X, 4), collapse=", "), "\n")
} else {
  cat("✗ Algorithm failed to produce results\n")
}

cat("\n=== FULL PIPELINE TEST COMPLETE ===\n")
