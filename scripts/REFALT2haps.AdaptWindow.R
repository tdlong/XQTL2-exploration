#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Adaptive Window Testing Script
# =============================================================================
# This script tests the adaptive window algorithm on a single genomic region
# It is separate from the production REFALT2haps pipeline
# Usage: Rscript scripts/REFALT2haps.AdaptWindow.R chr parfile mydir test_pos test_window_size

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.R chr parfile mydir test_pos test_window_size\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.R chr3R helpfiles/haplotype_parameters.R process/test 10000000 2000000\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
test_pos <- as.numeric(args[4])
test_window_size <- as.numeric(args[5])

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")
rdsfile <- paste0(mydir, "/df3.", mychr, ".RDS")
fileout <- paste0(mydir, "/R.haps.", mychr, ".adaptive.RDS")

cat("=== REFALT2haps Adaptive Window Testing ===\n")
cat("Chromosome:", mychr, "\n")
cat("Test position:", test_pos, "\n")
cat("Test window size:", test_window_size, "bp\n")
cat("Input file:", filein, "\n")
cat("Output file:", fileout, "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Build constraint matrix for lsei based on founder groups
#' @param founder_groups List of founder groups from clustering
#' @param founder_names Vector of founder names in order
#' @param previous_constraints List of previous constraint matrices (E, F)
#' @return List with E matrix and F vector for lsei
build_constraint_matrix <- function(founder_groups, founder_names, previous_constraints = NULL) {
  n_founders <- length(founder_names)
  
  # Start with sum-to-1 constraint
  E <- matrix(rep(1, n_founders), nrow = 1)
  F <- 1.0
  
  # Add constraints for each group
  for (group in founder_groups) {
    if (length(group) == 1) {
      # Individual founder - add equality constraint if we have previous estimate
      founder_idx <- which(founder_names == group[1])
      if (length(founder_idx) > 0) {
        constraint_row <- rep(0, n_founders)
        constraint_row[founder_idx] <- 1
        E <- rbind(E, constraint_row)
        F <- c(F, previous_constraints$individual_freqs[group[1]])
      }
    } else {
      # Group of founders - add sum constraint
      constraint_row <- rep(0, n_founders)
      for (founder in group) {
        founder_idx <- which(founder_names == founder)
        if (length(founder_idx) > 0) {
          constraint_row[founder_idx] <- 1
        }
      }
      E <- rbind(E, constraint_row)
      F <- c(F, previous_constraints$group_sums[paste(group, collapse = "+")])
    }
  }
  
  list(E = E, F = F)
}

#' Estimate haplotype for a single sample with adaptive constraints
#' @param sample_data SNP data for one sample
#' @param window_size Current window size being processed
#' @param previous_constraints Constraints from smaller windows
#' @param current_h_cutoff Current h_cutoff value being tested
#' @return List with Groups, Haps, Err, Names
estimate_single_sample_haplotype_adaptive <- function(sample_data, window_size, previous_constraints = NULL, current_h_cutoff = NULL) {
  # Extract founder matrix and sample frequencies
  founder_matrix <- sample_data %>% select(matches(founders))
  sample_freqs <- sample_data$freq
  
  # Filter for non-NA values
  valid_positions <- !is.na(sample_freqs)
  sample_freqs <- sample_freqs[valid_positions]
  founder_matrix <- founder_matrix[valid_positions, ]
  
  # Convert to matrix for clustering and optimization
  founder_matrix <- as.matrix(founder_matrix)
  
  # Use the current h_cutoff if provided, otherwise use default
  if (is.null(current_h_cutoff)) {
    adaptive_h_cutoff <- h_cutoff * 2  # More reasonable multiplier
  } else {
    adaptive_h_cutoff <- current_h_cutoff
  }
  
  # Cluster founders based on similarity
  founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = adaptive_h_cutoff)
  
  # Build constraint matrix if we have previous constraints
  if (!is.null(previous_constraints)) {
    founder_names <- colnames(founder_matrix)
    constraint_info <- build_constraint_matrix(founder_clusters, founder_names, previous_constraints)
    E <- constraint_info$E
    F <- constraint_info$F
  } else {
    # First round - just sum-to-1 constraint
    n_founders <- ncol(founder_matrix)
    E <- matrix(rep(1, n_founders), nrow = 1)
    F <- 1.0
  }
  
  # Solve constrained least squares for haplotype frequencies
  n_founders <- ncol(founder_matrix)
  constraints <- list(
    A = founder_matrix,
    B = sample_freqs,
    E = E,                           # Equality constraints
    F = F,                           # Right-hand side values
    G = diag(n_founders),            # Non-negative constraint
    H = matrix(rep(0.0003, n_founders)) # Small minimum value
  )
  
  # Solve using limSolve
  solution <- lsei(
    A = constraints$A,
    B = constraints$B,
    E = constraints$E,
    F = constraints$F,
    G = constraints$G,
    H = constraints$H,
    verbose = TRUE,
    fulloutput = TRUE
  )
  
  # Store constraint information for next round
  constraint_info <- list(
    window_size = window_size,
    founder_groups = founder_clusters,
    individual_freqs = setNames(solution$X, colnames(founder_matrix)),
    group_sums = sapply(unique(founder_clusters), function(group_id) {
      group_founders <- names(founder_clusters[founder_clusters == group_id])
      sum(solution$X[colnames(founder_matrix) %in% group_founders])
    })
  )
  
  list(
    Groups = founder_clusters,
    Haps = solution$X,
    Err = solution$cov,
    Names = names(solution$X),
    ConstraintInfo = constraint_info
  )
}

#' Run adaptive window haplotype estimation on single test window
#' @param snp_data Processed SNP dataset
#' @param chromosome Current chromosome being processed
#' @param test_pos Test position (center of window)
#' @param test_window_size Test window size
#' @return List with results from all window sizes
run_adaptive_window_test <- function(snp_data, chromosome, test_pos, test_window_size) {
  
  # Define h_cutoffs to test
  h_cutoffs_to_test <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)
  
  # Define window sizes to try (in base pairs)
  window_sizes <- c(10000, 25000, 50000, 100000, 200000, 500000, test_window_size)
  
  # Initialize results storage for all h_cutoffs
  all_h_cutoff_results <- list()
  
  cat("Testing adaptive window algorithm with different h_cutoffs\n")
  cat("H_cutoffs to test:", paste(h_cutoffs_to_test, collapse = ", "), "\n")
  cat("Window sizes to try:", paste(window_sizes, collapse = " → "), "bp\n")
  cat("Test window centered at position:", test_pos, "\n\n")
  
  # Test each h_cutoff
  for (hc_idx in seq_along(h_cutoffs_to_test)) {
    current_h_cutoff <- h_cutoffs_to_test[hc_idx]
    cat("\n=== Testing h_cutoff:", current_h_cutoff, "===\n")
    
    # Initialize results storage for this h_cutoff
    all_results <- list()
    previous_constraints <- NULL
    
    # Run adaptive window algorithm for this h_cutoff
    for (window_idx in seq_along(window_sizes)) {
      current_window_size <- window_sizes[window_idx]
      cat("\n--- Processing window size:", current_window_size, "bp (h_cutoff:", current_h_cutoff, ") ---\n")
      
      # Create single test window centered at test position
      spots <- data.frame(
        CHROM = chromosome,
        start = test_pos - current_window_size/2,
        end = test_pos + current_window_size/2
      )
      
      # Ensure window doesn't extend beyond chromosome boundaries
      min_pos <- min(snp_data$POS)
      max_pos <- max(snp_data$POS)
      if (spots$start < min_pos || spots$end > max_pos) {
        cat("Warning: Test window extends beyond chromosome boundaries. Skipping.\n")
        next
      }
      
      # Filter SNPs in this window for founders and samples of interest
      window_snps <- snp_data %>%
        filter(CHROM == spots$CHROM &
               POS > spots$start &
               POS < spots$end &
               (name %in% founders | name %in% names_in_bam)) %>%
        select(-c(CHROM, N)) %>%
        pivot_wider(names_from = name, values_from = freq) %>%
        pivot_longer(!c("POS", matches(founders)), 
                    names_to = "sample", values_to = "freq") %>%
        select(-POS)
      
      if (nrow(window_snps) == 0) {
        cat("Warning: No SNPs found in window. Skipping.\n")
        next
      }
      
      # Estimate haplotypes for each sample with current constraints
      cat("Estimating haplotypes with current constraints...\n")
      sample_haplotypes <- window_snps %>%
        group_by(sample) %>%
        nest() %>%
        mutate(
          haplotypes = map(data, estimate_single_sample_haplotype_adaptive, 
                          current_window_size, previous_constraints, current_h_cutoff)
        ) %>%
        select(-data) %>%
        unnest_wider(haplotypes)
      
      # Store results for this window size
      all_results[[as.character(current_window_size)]] <- sample_haplotypes
      
      # Extract constraint information for next round
      if (nrow(sample_haplotypes) > 0) {
        # Get constraint info from first sample (assuming all samples have similar structure)
        first_sample <- sample_haplotypes$ConstraintInfo[[1]]
        if (!is.null(first_sample)) {
          previous_constraints <- first_sample
          cat("✓ Constraints updated for next window size\n")
          
          # Show current founder groupings
          cat("Current founder groups:\n")
          for (group_id in unique(first_sample$founder_groups)) {
            group_founders <- names(first_sample$founder_groups[first_sample$founder_groups == group_id])
            cat("  Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
          }
        }
      }
      
      # Check if we've achieved complete founder separation
      if (!is.null(previous_constraints)) {
        n_groups <- length(unique(previous_constraints$founder_groups))
        cat("Current number of founder groups:", n_groups, "out of", length(founders), "founders\n")
        
        if (n_groups == length(founders)) {
          cat("✓ All founders are in separate groups. Adaptive algorithm complete for h_cutoff", current_h_cutoff, "\n")
          break  # Stop when we have 8 groups for 8 founders
        }
      }
    }
    
    # Store results for this h_cutoff
    all_h_cutoff_results[[as.character(current_h_cutoff)]] <- all_results
    
    cat("\n=== Completed h_cutoff", current_h_cutoff, "===\n")
    cat("Processed", length(all_results), "window sizes\n")
  }
  
  cat("\n=== Adaptive window estimation complete ===\n")
  cat("Tested", length(h_cutoffs_to_test), "h_cutoff values\n")
  
  return(all_h_cutoff_results)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

cat("Loading and preprocessing data...\n")

# Load and preprocess REF/ALT count data
cat("Loading REF/ALT data...\n")
df <- read.table(filein, header = TRUE)

# Transform REF/ALT counts to frequencies
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

rm(df)  # Free memory
cat("✓ Frequency data prepared\n")

# Identify high-quality SNPs for analysis
cat("Filtering for high-quality SNPs...\n")
good_snps <- df2 %>%
  filter(name %in% founders) %>%
  group_by(CHROM, POS) %>%
  summarize(
    zeros = sum(N == 0),                    # No missing data
    not_fixed = sum(N != 0 & freq > 0.03 & freq < 0.97),  # Not fixed in founders
    informative = (sum(freq) > 0.05 | sum(freq) < 0.95)   # Sufficiently variable
  ) %>%
  ungroup() %>%
  filter(zeros == 0 & not_fixed == 0 & informative == TRUE) %>%
  select(c(CHROM, POS))

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

rm(df2)  # Free memory
cat("✓ High-quality SNP dataset created\n")

# Save intermediate dataset
saveRDS(df3, file = rdsfile)
cat("✓ Intermediate data saved to RDS\n")

# Run adaptive window test
cat("\nRunning adaptive window test...\n")
adaptive_results <- run_adaptive_window_test(
  snp_data = df3,
  chromosome = mychr,
  test_pos = test_pos,
  test_window_size = test_window_size
)

# Save results
saveRDS(adaptive_results, file = fileout)
cat("✓ Adaptive window results saved to:", fileout, "\n")

# Show summary of results
cat("\n=== RESULTS SUMMARY ===\n")
for (h_cutoff_key in names(adaptive_results)) {
  h_cutoff_results <- adaptive_results[[h_cutoff_key]]
  for (window_size in names(h_cutoff_results)) {
    result <- h_cutoff_results[[window_size]]
    if (nrow(result) > 0) {
      first_sample <- result$ConstraintInfo[[1]]
      if (!is.null(first_sample)) {
        n_groups <- length(unique(first_sample$founder_groups))
        cat("h_cutoff", h_cutoff_key, "Window size", window_size, "bp:", n_groups, "founder groups\n")
      }
    }
  }
}

cat("\nREFALT2haps Adaptive Window test completed successfully!\n")
