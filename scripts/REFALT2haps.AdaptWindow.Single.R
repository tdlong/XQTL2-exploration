#!/usr/bin/env Rscript

# REFALT2haps Adaptive Window - Haplotype Frequency Estimation with Constraint Accumulation
# This script estimates haplotype frequencies using progressive window expansion
# with constraint accumulation from smaller windows
# 
# Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>
# Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 4

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 4\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
h_cutoff_parameter <- as.numeric(args[4])

# Source the parameter file
source(parfile)

# Use h_cutoff from command line (ignores parameter file value)
h_cutoff_used <- h_cutoff_parameter

# Display parameters
cat("Parameters loaded from", parfile, ":\n")
cat("  Founders (", length(founders), "):", paste(founders, collapse=", "), "\n")
cat("  Step:", step, "\n")
cat("  H_cutoff (from command line):", h_cutoff_used, "\n")
cat("  Samples to process (", length(names_in_bam), "):", paste(names_in_bam, collapse=", "), "\n")

# Load and transform data (same as test script)
filein <- file.path(mydir, paste0("RefAlt.", mychr, ".txt"))
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

# Transform to wide format and apply quality filter ONCE
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

# Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions and include sample data
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

# Get non-founder samples from parameter file
non_founder_samples <- names_in_bam

# Set up scanning positions based on step parameter
scan_positions <- seq(from = min(quality_filtered_positions), 
                     to = max(quality_filtered_positions), 
                     by = step)

# Initialize results list
results_list <- list()

# Process each position and sample combination  
for (test_pos in scan_positions) {
  for (sample_name in non_founder_samples) {
    
    # ADAPTIVE WINDOW ALGORITHM WITH CONSTRAINT ACCUMULATION
    # (Same as main test above, but multiple cases)
    
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
    best_result <- NULL
    best_n_groups <- 0
    best_window_size <- window_sizes[1]
    previous_n_groups <- 0
    
    # CONSTRAINT ACCUMULATION (core of adaptive algorithm)
    accumulated_constraints <- NULL
    accumulated_constraint_values <- NULL
    
    for (window_idx in seq_along(window_sizes)) {
      window_size <- window_sizes[window_idx]
      window_start <- max(1, test_pos - window_size/2)
      window_end <- test_pos + window_size/2
      
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
      
      # Hierarchical clustering
      founder_dist <- dist(t(founder_matrix_clean))
      hclust_result <- hclust(founder_dist, method = "complete")
      groups <- cutree(hclust_result, h = h_cutoff_used)
      n_groups <- length(unique(groups))
      
      # Check if clustering improved
      if (window_idx > 1 && n_groups <= previous_n_groups) {
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
      }
      
      # Run LSEI with constraints
      tryCatch({
        result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                E = E, F = F, 
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
        
        if (result$IsError == 0) {
          # LSEI successful - accumulate constraints for next window
          current_constraints <- NULL
          current_constraint_values <- NULL
          
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
              
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
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
          } else {
            accumulated_constraints <- NULL
            accumulated_constraint_values <- NULL
          }
          
          # Store the best result
          best_result <- result
          best_n_groups <- n_groups
          best_window_size <- window_size
          
          # Check if all founders are separated
          if (n_groups == length(founders)) {
            break  # Success! Stop expanding
          }
        }
      }, error = function(e) {
        # LSEI error, try larger window
      })
    }
    
    # Apply the correct rules for output
    if (!is.null(best_result)) {
      # LSEI was successful - ALWAYS return the frequency estimates
      founder_frequencies <- best_result$X
      
      # Check if founders are distinguishable to set trust level
      if (best_n_groups == length(founders)) {
        estimate_OK <- 1  # Founders distinguishable - we can TRUST the frequencies
      } else {
        estimate_OK <- 0  # Founders NOT distinguishable - we CANNOT trust the frequencies
      }
      
    } else {
      # Either insufficient SNPs OR LSEI failed/didn't converge
      # Both cases → NA for estimate_OK and frequencies
      founder_frequencies <- rep(NA, length(founders))
      estimate_OK <- NA
    }
    
    # CREATE EXACT SAME result_row STRUCTURE AS PRODUCTION
    result_row <- list(
      chr = mychr,
      pos = test_pos,
      sample = sample_name,
      h_cutoff = h_cutoff_used,
      final_window_size = best_window_size,
      n_snps = ifelse(is.null(best_result), 0, nrow(wide_data)),
      estimate_OK = estimate_OK
    )
    
    # Add founder frequencies as named columns (EXACTLY like production should do)
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- founder_frequencies[i]
    }
    
    results_list[[length(results_list) + 1]] <- result_row
  }
}

# Convert to data frame and save
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(mydir, "haplotype_results")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  output_file <- file.path(output_dir, paste0("adaptive_window_h", h_cutoff_used, "_results_", mychr, ".RDS"))
  saveRDS(results_df, output_file)
  
  cat("Adaptive window haplotype estimation complete.\n")
  cat("Processed", nrow(results_df), "position/sample combinations\n")
  cat("Results saved to:", output_file, "\n")
  
} else {
  cat("No results generated - check data and parameters\n")
  quit(status = 1)
}