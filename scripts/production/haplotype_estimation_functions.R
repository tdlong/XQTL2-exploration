#!/usr/bin/env Rscript

# Haplotype Estimation Functions
# Professional implementation with unified algorithm and configurable verbosity

library(tidyverse)
library(limSolve)

#' Estimate Haplotype Frequencies (Adaptive Window Only)
#' 
#' Adaptive window haplotype estimation with hierarchical clustering
#' 
#' @param pos Genomic position (integer)
#' @param sample_name Sample identifier (string)
#' @param df3 Processed data frame with POS, name, freq columns
#' @param founders Vector of founder names
#' @param h_cutoff Hierarchical clustering cutoff threshold
#' @param method Method (ignored - always uses adaptive window)
#' @param window_size_bp Window size in base pairs (ignored - always uses adaptive window)
#' @param chr Chromosome name (for output)
#' @param verbose Verbosity level (0=silent, 1=basic, 2=detailed, 3=full debug)
#' 
#' @return List with:
#'   - Groups: vector of group assignments from clustering (length = number of founders)
#'   - Haps: vector of founder haplotype frequencies (length = number of founders, sums to 1.0)
#'   - Err: error covariance matrix from lsei (n_founders x n_founders)
#'   - Names: vector of founder names (length = number of founders)
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  # DEBUG: Confirm this modified function is being called
  cat("DEBUG: Using MODIFIED estimate_haplotypes function from working_copy.R\n")
  
  method <- match.arg(method)
  
  if (verbose >= 1) {
    cat(sprintf("Processing pos: %s, sample: %s, method: %s\n", 
                format(pos, big.mark=","), sample_name, method))
  }
    # ADAPTIVE WINDOW METHOD
    if (verbose >= 2) {
      cat(sprintf("=== ADAPTIVE WINDOW METHOD: h_cutoff = %g ===\n", h_cutoff))
      cat(sprintf("Position: %s\n", format(pos, big.mark=",")))
      cat(sprintf("Sample: %s\n", sample_name))
    }
    
    window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
    if (verbose >= 2) {
      cat("Window sizes to test:", paste(window_sizes/1000, "kb"), "\n")
    }
    
    final_result <- NULL
    final_n_groups <- 0
    final_window_size <- window_sizes[1]
    final_wide_data <- NULL
    previous_n_groups <- 0
    
    # CONSTRAINT ACCUMULATION (core of adaptive algorithm)
    accumulated_constraints <- NULL
    accumulated_constraint_values <- NULL
    
    for (window_idx in seq_along(window_sizes)) {
      window_size <- window_sizes[window_idx]
      window_start <- max(1, pos - window_size/2)
      window_end <- pos + window_size/2
      
      if (verbose >= 2) {
        cat(sprintf("\n--- Window %d: %s ---\n", window_idx, 
                    ifelse(window_size >= 1000, paste0(window_size/1000, " kb"), paste0(window_size, " bp"))))
      }
      
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
      
      # Always update final window info (this is the last window we actually tried)
      final_window_size <- window_size
      final_wide_data <- wide_data
      
      # Hierarchical clustering
      founder_dist <- dist(t(founder_matrix_clean))
      hclust_result <- hclust(founder_dist, method = "complete")
      groups <- cutree(hclust_result, h = h_cutoff)
      n_groups <- length(unique(groups))
      
      if (verbose >= 2) {
        cat(sprintf("SNPs in window: %d positions\n", nrow(founder_matrix_clean)))
        cat(sprintf("Clustering results: %d groups\n", n_groups))
        
        if (verbose >= 3 && n_groups <= 6) {
          # Show group assignments
          unique_clusters <- unique(groups)
          for (cluster_id in unique_clusters) {
            cluster_founders <- founders[groups == cluster_id]
            if (length(cluster_founders) == 1) {
              cat(sprintf("  Group %d: %s (individual)\n", cluster_id, cluster_founders))
            } else {
              cat(sprintf("  Group %d: %s\n", cluster_id, paste(cluster_founders, collapse=", ")))
            }
          }
        }
      }
      
      # Check if clustering improved
      if (window_idx > 1 && n_groups <= previous_n_groups) {
        if (verbose >= 2) cat("✗ No improvement from previous window, trying larger window\n")
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
        if (verbose >= 2) {
          cat(sprintf("    ✓ Added %d accumulated constraints from previous windows\n", nrow(accumulated_constraints)))
        }
      } else {
        if (verbose >= 2) {
          cat("    • No accumulated constraints yet (first meaningful window)\n")
        }
      }
      
      # Run LSEI with constraints
      tryCatch({
        result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                E = E, F = F, 
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE)
        
        if (result$IsError == 0) {
          # LSEI successful - accumulate constraints for next window
          current_constraints <- NULL
          current_constraint_values <- NULL
          
          if (verbose >= 2) {
            cat(sprintf("    Building constraints from %d founder groups...\n", n_groups))
          }
          
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
              
              if (verbose >= 2) {
                group_names <- founders[cluster_founders]
                cat(sprintf("      Group %d: %s = %.4f\n", 
                           cluster_id, paste(group_names, collapse="+"), group_freq))
              }
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
              # Create constraint: this founder = their exact frequency
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
              
              if (verbose >= 2) {
                cat(sprintf("      Group %d: %s = %.4f (locked)\n", 
                           cluster_id, founders[cluster_founders], founder_freq))
              }
            }
          }
          
          # Update accumulated constraints for next window
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
            if (verbose >= 2) {
              cat(sprintf("    ✓ Accumulated %d constraints for next window (groups: %s)\n", 
                         nrow(current_constraints),
                         paste(unique(groups), collapse=",")))
            }
          } else {
            accumulated_constraints <- NULL
            accumulated_constraint_values <- NULL
            if (verbose >= 2) {
              cat("    • No constraints to accumulate (all founders in 1 group)\n")
            }
          }
          
          # Store the result (will be overwritten as we expand)
          final_result <- result
          final_n_groups <- n_groups
          
          # Check if all founders are separated
          if (n_groups == length(founders)) {
            if (verbose >= 2) {
              cat("✅ SUCCESS: All founders distinguished - stopping window expansion\n")
            }
            break  # Success! Stop expanding
          }
        }
      }, error = function(e) {
        # LSEI error, try larger window
        if (verbose >= 2) {
          cat(sprintf("LSEI error: %s\n", e$message))
        }
      })
    }
    
    # Apply the correct rules for output
    if (!is.null(final_result)) {
      # LSEI was successful - ALWAYS return the frequency estimates
      founder_frequencies <- final_result$X
      names(founder_frequencies) <- founders
      
      # Check if founders are distinguishable to set trust level
      if (final_n_groups == length(founders)) {
        estimate_OK <- 1  # Founders distinguishable - we can TRUST the frequencies
        if (verbose >= 2) {
          cat("✓ Founders distinguishable - frequencies are TRUSTWORTHY (estimate_OK = 1)\n")
        }
      } else {
        estimate_OK <- 0  # Founders NOT distinguishable - we CANNOT trust the frequencies
        if (verbose >= 2) {
          cat("⚠️  Founders NOT fully distinguished - frequencies are UNTRUSTWORTHY (estimate_OK = 0)\n")
        }
      }
      
    } else {
      # Either insufficient SNPs OR LSEI failed/didn't converge
      founder_frequencies <- rep(NA, length(founders))
      names(founder_frequencies) <- founders
      estimate_OK <- NA
      if (verbose >= 1) {
        cat("❌ LSEI failed or insufficient SNPs - no haplotype frequencies available\n")
      }
    }
    
    if (verbose >= 2) {
      cat(sprintf("\n✅ FINAL RESULT:\n"))
      cat(sprintf("  estimate_OK: %s\n", ifelse(is.na(estimate_OK), "NA", estimate_OK)))
      cat(sprintf("  final_window_size: %s\n", 
                  ifelse(final_window_size >= 1000, paste0(final_window_size/1000, " kb"), paste0(final_window_size, " bp"))))
      cat(sprintf("  n_snps: %d\n", ifelse(is.null(final_result), 0, nrow(final_wide_data))))
    }
    
    return(list(Groups=groups, Haps=founder_frequencies, Err=result$cov, Names=founders))
  
}


#' Show SNP diagnostic information
show_snp_diagnostics <- function(wide_data, founders, founder_matrix_clean) {
  cat("\n=== FOUNDER FREQUENCY RANGES ===\n")
  for (i in seq_along(founders)) {
    founder_name <- founders[i]
    freq_range <- range(founder_matrix_clean[, i], na.rm = TRUE)
    cat(founder_name, ": ", sprintf("%.3f - %.3f", freq_range[1], freq_range[2]), "\n")
  }
  
  cat("\n=== SAMPLE SNP DATA (first 10 positions, percentages) ===\n")
  sample_data <- wide_data[1:min(10, nrow(wide_data)), ]
  
  # Create formatted table
  cat(sprintf("%-10s", "POS"), paste(sprintf("%-3s", founders), collapse=" "), "\n")
  cat(paste(rep("-", 10 + length(founders) * 4), collapse=""), "\n")
  
  for (i in 1:nrow(sample_data)) {
    pos <- sample_data$POS[i]
    freqs <- sprintf("%2.0f", as.numeric(sample_data[i, founders]) * 100)
    cat(sprintf("%-10s", pos), paste(sprintf("%3s", freqs), collapse=" "), "\n")
  }
}

#' Show clustering diagnostic information  
show_clustering_diagnostics <- function(distances, founders, h_cutoff, groups, n_groups, haplotype_freqs) {
  cat("\n=== CLUSTERING ANALYSIS ===\n")
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
  
  cat("\nClustering results at h_cutoff", h_cutoff, ":\n")
  cat("Number of groups:", n_groups, "\n")
  cat("Group assignments:", paste(groups, collapse=", "), "\n")
  
  # Show which founders are in which groups
  for (group_id in unique(groups)) {
    group_founders <- founders[groups == group_id]
    cat("Group", group_id, ":", paste(group_founders, collapse=", "), "\n")
  }
  
  # Show group-wise frequencies
  cat("\nGroup-wise frequencies:\n")
  for (group_id in unique(groups)) {
    group_founders <- founders[groups == group_id]
    group_freq <- sum(haplotype_freqs[groups == group_id])
    cat("Group", group_id, ":", paste(group_founders, collapse="+"), "=", round(group_freq, 4), "\n")
  }
}
