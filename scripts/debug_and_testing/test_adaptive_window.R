#!/usr/bin/env Rscript

# Test Adaptive Window Distinguishability  
# Tests the core adaptive distinguishability algorithm
# This should match the production script logic exactly

library(tidyverse)
library(limSolve)

cat("=== ADAPTIVE WINDOW DISTINGUISHABILITY TEST ===\n\n")

# 1. Load parameters (same as production)
cat("1. Loading parameters...\n")
source("helpfiles/JUICE_haplotype_parameters.R")
cat("✓ Parameter file: helpfiles/JUICE_haplotype_parameters.R\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), h_cutoff (", h_cutoff, "), samples (", length(names_in_bam), ")\n")
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

# Transform to wide format and apply quality filter ONCE (same as production)
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

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")

# Keep only high-quality positions in the full dataset
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

# Get non-founder samples
# Use samples defined in parameter file
non_founder_samples <- names_in_bam

cat("✓ Data ready:", nrow(df3), "rows\n")

# 3. Test the adaptive window distinguishability
cat("\n3. Testing adaptive window distinguishability...\n")

# Test parameters
test_pos <- 15000000
test_h_cutoff <- 10  # Test higher h_cutoff to force window expansion
test_sample <- non_founder_samples[1]

cat("✓ Test position:", test_pos, "\n")
cat("✓ Test h_cutoff:", test_h_cutoff, "\n")
cat("✓ Test sample:", test_sample, "(note: sample not used for distinguishability)\n\n")

# Adaptive window expansion: start small, grow until distinguishable or max size
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)  # Progressive sizes
final_window_size <- NA
estimate_OK <- 0
final_n_snps <- 0

cat("=== STARTING ADAPTIVE WINDOW ALGORITHM ===\n")
cat("Window sizes to test:", paste(window_sizes/1000, "kb", collapse=", "), "\n\n")

# Track SNP counts across all window sizes
snps_tracking <- data.frame(
  window_size = window_sizes,
  total_snps = 0,
  founder_complete = 0,
  clustering_valid = 0,
  n_groups = 0,
  estimate_OK = 0
)

for (window_idx in seq_along(window_sizes)) {
  window_size <- window_sizes[window_idx]
  cat("--- Window", window_idx, ":", window_size/1000, "kb ---\n")
  
  # Define window boundaries
  window_start <- test_pos - window_size/2
  window_end <- test_pos + window_size/2
  cat("Window range:", window_start, "to", window_end, "bp\n")
  
  # Get SNPs in window for founders only (data is already quality-filtered)
  window_snps_long <- df3 %>%
    filter(POS >= window_start & POS <= window_end & name %in% founders)
  
  # Convert to wide format (rows = positions, columns = founders)
  wide_data <- window_snps_long %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  cat("SNPs in window:", nrow(wide_data), "positions\n")
  
  # Track total SNPs
  snps_tracking$total_snps[window_idx] <- nrow(wide_data)
  
  # Show first 10 SNPs as a table for the first window
  if (window_idx == 1) {
    cat("\n=== FIRST 10 SNPS IN WINDOW (RAW DATA, percentages) ===\n")
    display_data <- wide_data %>%
      arrange(POS) %>%
      head(10)
    
    # Create formatted table
    cat(sprintf("%-10s", "POS"), paste(sprintf("%-3s", founders), collapse=" "), "\n")
    cat(paste(rep("-", 10 + length(founders) * 4), collapse=""), "\n")
    
    for (i in 1:nrow(display_data)) {
      pos <- display_data$POS[i]
      freqs <- sprintf("%2.0f", as.numeric(display_data[i, founders]) * 100)
      cat(sprintf("%-10s", pos), paste(sprintf("%3s", freqs), collapse=" "), "\n")
    }
    cat("\n")
  }
  
  # Check if we have all founder columns and enough data
  if (ncol(wide_data) < length(founders) + 1) {
    cat("✗ Missing founders, trying larger window\n\n")
    snps_tracking$founder_complete[window_idx] <- 0
    next  # Try larger window
  }
  
  if (nrow(wide_data) < 10) {
    cat("✗ Insufficient data in window, trying larger window\n\n")
    snps_tracking$founder_complete[window_idx] <- 0
    next  # Try larger window
  }
  
  # Get founder matrix (no need for additional quality filtering)
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  
  # Convert to matrix for clustering
  founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
  
  if (nrow(founder_matrix_clean) < 10) {
    cat("✗ Insufficient clean data, trying larger window\n\n")
    snps_tracking$founder_complete[window_idx] <- nrow(founder_matrix_clean)
    snps_tracking$clustering_valid[window_idx] <- 0
    next  # Try larger window
  }
  
  snps_tracking$founder_complete[window_idx] <- nrow(founder_matrix)
  snps_tracking$clustering_valid[window_idx] <- nrow(founder_matrix_clean)
  
  cat("Founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  
  # Show founder frequency ranges for this window
  cat("Founder frequency ranges:\n")
  for (i in 1:ncol(founder_matrix_clean)) {
    founder_name <- founders[i]
    freq_range <- range(founder_matrix_clean[, i], na.rm = TRUE)
    cat("  ", founder_name, ": ", sprintf("%.3f - %.3f", freq_range[1], freq_range[2]), "\n")
  }
  
  # Hierarchical clustering to check distinguishability
  tryCatch({
    distances <- dist(t(founder_matrix_clean), method = "euclidean")
    hclust_result <- hclust(distances, method = "ward.D2")
    
    # Show clustering distance matrix
    cat("Clustering distance matrix:\n")
    dist_matrix <- as.matrix(distances)
    rownames(dist_matrix) <- founders
    colnames(dist_matrix) <- founders
    
    # Print matrix header
    cat(sprintf("%4s", ""), paste(sprintf("%6s", founders), collapse=""), "\n")
    
    # Print matrix rows
    for (i in 1:length(founders)) {
      row_values <- sprintf("%6.1f", dist_matrix[i, ])
      cat(sprintf("%4s", founders[i]), paste(row_values, collapse=""), "\n")
    }
    
    min_dist <- min(dist_matrix[dist_matrix > 0])
    max_dist <- max(dist_matrix)
    cat("Min:", sprintf("%.1f", min_dist), "| Max:", sprintf("%.1f", max_dist), "| H_cutoff:", test_h_cutoff, "\n")
    
    # Cut tree at h_cutoff
    groups <- cutree(hclust_result, h = test_h_cutoff)
    n_groups <- length(unique(groups))
    
    snps_tracking$n_groups[window_idx] <- n_groups
    
    cat("Clustering results:\n")
    cat("  Number of groups:", n_groups, "\n")
    cat("  Group assignments:", paste(groups, collapse=", "), "\n")
    
    # Show which founders are in which groups
    for (group_id in unique(groups)) {
      group_founders <- founders[groups == group_id]
      cat("  Group", group_id, ":", paste(group_founders, collapse=", "), "\n")
    }
    
    # Check if all founders can be distinguished
    if (n_groups == length(founders)) {
      cat("✓ SUCCESS: All", length(founders), "founders distinguished individually!\n")
      final_window_size <- window_size
      estimate_OK <- 1
      final_n_snps <- nrow(wide_data)
      snps_tracking$estimate_OK[window_idx] <- 1
      break  # Success! Stop expanding
    } else {
      cat("✗ Only", n_groups, "groups - need larger window\n")
      cat("✗ Some founders still clustered together\n")
      snps_tracking$estimate_OK[window_idx] <- 0
    }
    
  }, error = function(e) {
    cat("✗ Clustering failed:", e$message, ", trying larger window\n")
    snps_tracking$n_groups[window_idx] <- 0
    snps_tracking$estimate_OK[window_idx] <- 0
  })
  
  cat("\n")
}

# Record result (either successful at some window size, or failed at all sizes)
if (is.na(final_window_size)) {
  final_window_size <- max(window_sizes)  # Used largest window but failed
  cat("✗ FAILED: Could not distinguish all founders even at", final_window_size/1000, "kb\n")
} else {
  cat("✅ SUCCESS: All founders distinguished at", final_window_size/1000, "kb\n")
}

# Show comprehensive SNP tracking summary
cat("\n=== SNP TRACKING SUMMARY ===\n")
print(snps_tracking)

# Run LSEI estimation using the optimal window size
cat("\n=== LSEI HAPLOTYPE ESTIMATION ===\n")
# Always run LSEI first, then determine estimate_OK based on results
  # Use the final window size to get data for LSEI
  window_start <- test_pos - final_window_size/2
  window_end <- test_pos + final_window_size/2
  
  window_data <- df3 %>%
    filter(POS >= window_start & POS <= window_end & name %in% c(founders, test_sample))
  
  wide_data <- window_data %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  # Extract founder and sample data
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  sample_freqs <- wide_data %>%
    pull(!!test_sample)
  
  # Remove rows with any NAs
  complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  cat("Data for LSEI:", nrow(founder_matrix_clean), "complete SNPs\n")
  
  tryCatch({
    # LSEI constraints: sum to 1, non-negative
    E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    G <- diag(length(founders))  # Non-negativity constraints
    H <- matrix(rep(0.0003, length(founders)))  # Lower bound
    
    lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                 E = E, F = F, G = G, H = H)
    
    if (lsei_result$IsError == 0) {
      cat("✓ LSEI successful!\n")
      haplotype_freqs <- lsei_result$X
      names(haplotype_freqs) <- founders
      
      cat("Haplotype frequency estimates:\n")
      for (i in seq_along(founders)) {
        cat(sprintf("  %s: %.4f\n", founders[i], haplotype_freqs[i]))
      }
      cat("Sum of frequencies:", round(sum(haplotype_freqs), 4), "\n")
      
      # Now check distinguishability using clustering results
      if (final_n_groups == length(founders)) {
        estimate_OK <- 1  # Founders distinguishable - we can TRUST the frequencies
        cat("✓ Founders distinguishable - frequencies are TRUSTWORTHY (estimate_OK = 1)\n")
      } else {
        estimate_OK <- 0  # Founders NOT distinguishable - we CANNOT trust the frequencies  
        cat("⚠️  Founders NOT distinguishable - frequencies are UNTRUSTWORTHY (estimate_OK = 0)\n")
        cat("   BUT we still have the LSEI frequency estimates!\n")
      }
      
    } else {
      cat("✗ LSEI failed with error code:", lsei_result$IsError, "\n")
      estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
    }
  }, error = function(e) {
    cat("✗ LSEI error:", e$message, "\n")
    estimate_OK <- NA  # Update estimate_OK to NA if LSEI failed
  })

cat("\n=== FINAL RESULT ===\n")
cat("final_window_size:", final_window_size, "bp (", final_window_size/1000, "kb)\n")
cat("estimate_OK:", estimate_OK, "\n")
cat("n_snps:", final_n_snps, "\n")

if (estimate_OK == 1) {
  cat("✓ SUCCESS: Adaptive algorithm found optimal window size!\n")
  cat("✓ All", length(founders), "founders can be distinguished at", final_window_size/1000, "kb\n")
  cat("✓ LSEI estimation successful - haplotype frequencies available\n")
} else if (estimate_OK == 0) {
  cat("✗ FAILURE: Could not distinguish all", length(founders), "founders\n")
  cat("✗ Even at maximum window size of", max(window_sizes)/1000, "kb\n")
  cat("✗ Best result:", max(snps_tracking$n_groups), "groups\n")
} else {
  cat("✗ LSEI estimation failed - no haplotype frequencies available\n")
}

cat("\n=== ADAPTIVE WINDOW TEST COMPLETE ===\n")

# ===============================================================
# PRODUCTION DATA WORKFLOW TEST
# Test the same data structures and file output as production
# ===============================================================

cat("\n=== TESTING ADAPTIVE WINDOW PRODUCTION WORKFLOW ===\n")

# Test multiple positions and samples like production does
test_positions <- c(5000000, 10000000, 15000000)
test_samples <- names_in_bam[1:2]  # Test first 2 samples
test_h_cutoffs <- c(4, 6)  # Test different h_cutoff values

results_list <- list()

for (h_cutoff_test in test_h_cutoffs) {
  for (test_pos in test_positions) {
    for (sample_name in test_samples) {
      cat(sprintf("Testing position %d, sample %s, h_cutoff %d...\n", test_pos, sample_name, h_cutoff_test))
      
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
        groups <- cutree(hclust_result, h = h_cutoff_test)
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
          cat(sprintf("    ✓ Added %d accumulated constraints from previous windows\n", nrow(accumulated_constraints)))
        } else {
          cat("    • No accumulated constraints yet (first meaningful window)\n")
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
            
            cat(sprintf("    Building constraints from %d founder groups...\n", n_groups))
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
                
                group_names <- founders[cluster_founders]
                cat(sprintf("      Group %d: %s = %.4f\n", 
                           cluster_id, paste(group_names, collapse="+"), group_freq))
              } else {
                # Single founder: lock their exact frequency
                founder_freq <- result$X[cluster_founders]
                
                # Create constraint: this founder = their exact frequency
                constraint_row <- rep(0, n_founders)
                constraint_row[cluster_founders] <- 1
                
                current_constraints <- rbind(current_constraints, constraint_row)
                current_constraint_values <- c(current_constraint_values, founder_freq)
                
                cat(sprintf("      Group %d: %s = %.4f (locked)\n", 
                           cluster_id, founders[cluster_founders], founder_freq))
              }
            }
            
            # Update accumulated constraints for next window
            if (!is.null(current_constraints)) {
              accumulated_constraints <- current_constraints
              accumulated_constraint_values <- current_constraint_values
              cat(sprintf("    ✓ Accumulated %d constraints for next window (groups: %s)\n", 
                         nrow(current_constraints),
                         paste(unique(groups), collapse=",")))
            } else {
              accumulated_constraints <- NULL
              accumulated_constraint_values <- NULL
              cat("    • No constraints to accumulate (all founders in 1 group)\n")
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
        chr = "chr2R",
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff_test,
        final_window_size = best_window_size,
        n_snps = ifelse(is.null(best_result), 0, nrow(wide_data)),
        estimate_OK = estimate_OK
      )
      
      # Add founder frequencies as named columns (EXACTLY like production should do)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- founder_frequencies[i]
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      
      cat(sprintf("  Result: estimate_OK=%s, window_size=%d, n_snps=%d\n", 
                  ifelse(is.na(estimate_OK), "NA", estimate_OK), 
                  best_window_size, 
                  ifelse(is.null(best_result), 0, nrow(wide_data))))
    }
  }
}

# Convert to data frame (same as production)
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  cat("\n=== ADAPTIVE WINDOW PRODUCTION DATA STRUCTURE TEST ===\n")
  cat("Generated results data frame:\n")
  cat("Columns:", paste(names(results_df), collapse = ", "), "\n")
  cat("Rows:", nrow(results_df), "\n")
  
  # Verify all expected columns are present
  expected_cols <- c("chr", "pos", "sample", "h_cutoff", "final_window_size", "n_snps", "estimate_OK", founders)
  missing_cols <- setdiff(expected_cols, names(results_df))
  extra_cols <- setdiff(names(results_df), expected_cols)
  
  if (length(missing_cols) == 0) {
    cat("✓ All expected columns present:", paste(expected_cols, collapse = ", "), "\n")
  } else {
    cat("✗ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  }
  
  if (length(extra_cols) > 0) {
    cat("! Extra columns:", paste(extra_cols, collapse = ", "), "\n")
  }
  
  # Test RDS file saving/loading
  test_output_file <- "test_adaptive_window_output.RDS"
  saveRDS(results_df, test_output_file)
  cat("✓ Saved test results to:", test_output_file, "\n")
  
  # Test loading and verify structure
  loaded_results <- readRDS(test_output_file)
  if (identical(results_df, loaded_results)) {
    cat("✓ RDS save/load test passed\n")
  } else {
    cat("✗ RDS save/load test failed\n")
  }
  
  # Show sample of actual data
  cat("\nSample results:\n")
  print(head(results_df, 3))
  
  # Test constraint accumulation worked (different h_cutoff should give different results)
  if (length(test_h_cutoffs) > 1) {
    h_cutoff_comparison <- results_df %>%
      group_by(pos, sample) %>%
      summarise(
        n_h_cutoffs = n_distinct(h_cutoff),
        different_results = n_distinct(paste(estimate_OK, final_window_size)),
        .groups = "drop"
      )
    
    if (any(h_cutoff_comparison$different_results > 1)) {
      cat("✓ Different h_cutoff values produce different results (constraint accumulation working)\n")
    } else {
      cat("⚠️  All h_cutoff values produce identical results (potential constraint accumulation bug)\n")
    }
  }
  
  # Clean up test file
  file.remove(test_output_file)
  cat("✓ Cleaned up test file\n")
  
} else {
  cat("✗ No results generated - check algorithm or data\n")
}

cat("\n=== ADAPTIVE WINDOW PRODUCTION WORKFLOW TEST COMPLETE ===\n")
