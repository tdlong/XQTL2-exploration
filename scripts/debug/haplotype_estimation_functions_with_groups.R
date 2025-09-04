#!/usr/bin/env Rscript

# Haplotype Estimation Functions
# Professional implementation with unified algorithm and configurable verbosity

library(tidyverse)
library(limSolve)

#' Estimate Haplotype Frequencies
#' 
#' Unified function for fixed window and adaptive window haplotype estimation
#' 
#' @param pos Genomic position (integer)
#' @param sample_name Sample identifier (string)
#' @param df3 Processed data frame with POS, name, freq columns
#' @param founders Vector of founder names
#' @param h_cutoff Hierarchical clustering cutoff threshold
#' @param method Either "fixed" or "adaptive"
#' @param window_size_bp Window size in base pairs (required for fixed method)
#' @param chr Chromosome name (for output)
#' @param verbose Verbosity level (0=silent, 1=basic, 2=detailed, 3=full debug)
#' 
#' @return List with chr, pos, sample, method, final_window_size, n_snps, 
#'         estimate_OK, and individual founder frequencies
estimate_haplotypes <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = c("fixed", "adaptive"),
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  method <- match.arg(method)
  
  if (verbose >= 1) {
    cat(sprintf("Processing pos: %s, sample: %s, method: %s\n", 
                format(pos, big.mark=","), sample_name, method))
  }
  
  # Validate inputs
  if (method == "fixed" && is.null(window_size_bp)) {
    stop("window_size_bp required for fixed method")
  }
  
  if (method == "fixed") {
    # FIXED WINDOW METHOD
    if (verbose >= 2) {
      cat(sprintf("=== FIXED WINDOW METHOD: %d bp ===\n", window_size_bp))
      cat(sprintf("Position: %s\n", format(pos, big.mark=",")))
      cat(sprintf("Sample: %s\n", sample_name))
    }
    
    # Calculate window boundaries
    window_start <- pos - window_size_bp/2
    window_end <- pos + window_size_bp/2
    
    if (verbose >= 2) {
      cat(sprintf("Window range: %s to %s bp\n", 
                  format(window_start, big.mark=","), 
                  format(window_end, big.mark=",")))
    }
    
    # Get data (same as test script)
    window_data <- df3 %>%
      filter(POS >= window_start & POS <= window_end & name %in% c(founders, sample_name))
    
    if (nrow(window_data) == 0) {
      if (verbose >= 1) cat("No data in window\n")
      return(create_empty_result(chr, pos, sample_name, method, window_size_bp, founders))
    }
    
    wide_data <- window_data %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) {
      if (verbose >= 1) cat("Insufficient data in window\n")
      return(create_empty_result(chr, pos, sample_name, method, window_size_bp, founders))
    }
    
    if (verbose >= 2) {
      cat(sprintf("SNPs in window: %d positions\n", nrow(wide_data)))
    }
    
    # Get founder matrix and sample frequencies
    founder_matrix <- wide_data %>%
      select(all_of(founders)) %>%
      as.matrix()
    sample_freqs <- wide_data %>%
      pull(!!sample_name)
    
    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) {
      if (verbose >= 1) cat("Too few complete SNPs\n")
      return(create_empty_result(chr, pos, sample_name, method, window_size_bp, founders))
    }
    
    if (verbose >= 2) {
      cat(sprintf("Complete SNPs for analysis: %d\n", nrow(founder_matrix_clean)))
    }
    
    # Show diagnostic data if requested
    if (verbose >= 3) {
      show_snp_diagnostics(wide_data, founders, founder_matrix_clean)
    }
    
    # Run LSEI and clustering
    result <- run_lsei_and_clustering(founder_matrix_clean, sample_freqs_clean, 
                                     founders, h_cutoff, verbose)
    
    # Return result
    return(create_result(chr, pos, sample_name, method, window_size_bp, 
                        nrow(wide_data), result$estimate_OK, result$haplotype_freqs, founders,
                        groups = result$groups, error_matrix = result$error_matrix))
    
  } else {
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
    final_groups <- NULL
    final_error_matrix <- NULL
    
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
                                G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                                fulloutput = TRUE)  # Capture error matrix
        
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
          
          # Capture groups and error matrix
          final_groups <- groups
          names(final_groups) <- founders
          if (!is.null(result$cov)) {
            final_error_matrix <- result$cov
            rownames(final_error_matrix) <- founders
            colnames(final_error_matrix) <- founders
          }
          
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
    
    # Return result
    return(create_result(chr, pos, sample_name, method, final_window_size, 
                        ifelse(is.null(final_result), 0, nrow(final_wide_data)), 
                        estimate_OK, founder_frequencies, founders, h_cutoff,
                        groups = final_groups, error_matrix = final_error_matrix))
  }
}

#' Run LSEI estimation and clustering analysis
#' @param founder_matrix_clean Clean founder matrix (no NAs)
#' @param sample_freqs_clean Clean sample frequencies (no NAs) 
#' @param founders Vector of founder names
#' @param h_cutoff Clustering threshold
#' @param verbose Verbosity level
#' @return List with estimate_OK and haplotype_freqs
run_lsei_and_clustering <- function(founder_matrix_clean, sample_freqs_clean, 
                                   founders, h_cutoff, verbose) {
  
  # Initialize result variables
  estimate_OK <- NA
  haplotype_freqs <- rep(NA, length(founders))
  names(haplotype_freqs) <- founders
  
  # Initialize groups and error matrix variables
  groups <- rep(1, length(founders))
  names(groups) <- founders
  error_matrix <- matrix(NA, length(founders), length(founders))
  rownames(error_matrix) <- founders
  colnames(error_matrix) <- founders
  
  # ADDED: Detailed debugging output right before LSEI
  if (verbose >= 3) {
    cat("\n=== DETAILED DEBUGGING BEFORE LSEI ===\n")
    cat("Founder matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
    cat("Sample frequencies length:", length(sample_freqs_clean), "\n")
    cat("h_cutoff:", h_cutoff, "\n")
    
    # Show founder matrix summary
    cat("\nFounder matrix summary:\n")
    for (i in seq_along(founders)) {
      founder_data <- founder_matrix_clean[, i]
      cat(sprintf("  %s: %d SNPs, range %.3f-%.3f, mean %.3f\n", 
                 founders[i], length(founder_data), 
                 min(founder_data, na.rm=TRUE), max(founder_data, na.rm=TRUE),
                 mean(founder_data, na.rm=TRUE)))
    }
    
    # Show sample frequency summary
    cat("\nSample frequency summary:\n")
    cat("  Range:", min(sample_freqs_clean, na.rm=TRUE), "-", max(sample_freqs_clean, na.rm=TRUE), "\n")
    cat("  Mean:", mean(sample_freqs_clean, na.rm=TRUE), "\n")
    cat("  NAs:", sum(is.na(sample_freqs_clean)), "\n")
    
    # Calculate and show founder distances
    cat("\nFounder distance matrix:\n")
    founder_dist <- dist(t(founder_matrix_clean))
    founder_dist_matrix <- as.matrix(founder_dist)
    print(round(founder_dist_matrix, 4))
    
    # Show clustering results
    cat("\nHierarchical clustering results:\n")
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    
    cat("  h_cutoff:", h_cutoff, "\n")
    cat("  Number of groups:", n_groups, "\n")
    cat("  Expected groups:", length(founders), "\n")
    cat("  Groups sufficient:", n_groups == length(founders), "\n")
    
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
    
    # Show minimum distance between any two founders
    min_dist <- min(founder_dist)
    cat("\nDistance analysis:\n")
    cat("  Minimum distance between any two founders:", round(min_dist, 4), "\n")
    cat("  h_cutoff satisfied:", min_dist >= h_cutoff, "\n")
    
    # Show which founders are closest
    min_dist_idx <- which(founder_dist_matrix == min_dist, arr.ind = TRUE)
    if (nrow(min_dist_idx) > 0) {
      closest_pair <- min_dist_idx[1, ]
      cat("  Closest founders:", founders[closest_pair[1]], "and", founders[closest_pair[2]], "\n")
    }
    
    cat("=== END DEBUGGING ===\n\n")
  }
  
  # Run LSEI with better error handling
  tryCatch({
    E <- matrix(rep(1, length(founders)), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    G <- diag(length(founders))  # Non-negativity constraints
    H <- matrix(rep(0.0003, length(founders)))  # Lower bound
    
    if (verbose >= 2) {
      cat("\nLSEI inputs:\n")
      cat("Matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
      cat("Number of constraints:", nrow(E), "equality +", nrow(G), "inequality\n")
    }
    
    # Try LSEI with warning handler
    withCallingHandlers({
      lsei_result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                                   E = E, F = F, G = G, H = H)
    }, warning = function(w) {
      if (verbose >= 1) {
        cat("LSEI Warning:", conditionMessage(w), "\n")
      }
    })
    
    if (lsei_result$IsError == 0) {
      # LSEI successful - get frequencies
      haplotype_freqs <- lsei_result$X
      names(haplotype_freqs) <- founders
      
      if (verbose >= 1) {
        cat("✓ LSEI successful\n")
      }
      
      if (verbose >= 3) {
        cat("✓ LSEI successful!\n")
        cat("Haplotype frequency estimates:\n")
        for (i in seq_along(founders)) {
          cat(sprintf("  %s: %.4f\n", founders[i], haplotype_freqs[i]))
        }
        cat("Sum of frequencies:", round(sum(haplotype_freqs), 4), "\n")
      }
      
      # Check distinguishability
      distances <- dist(t(founder_matrix_clean))
      hclust_result <- hclust(distances, method = "complete")
      groups <- cutree(hclust_result, h = h_cutoff)
      names(groups) <- founders  # Add founder names to groups
      n_groups <- length(unique(groups))
      
      # Capture error matrix from LSEI result
      if (!is.null(lsei_result$cov)) {
        error_matrix <- lsei_result$cov
        rownames(error_matrix) <- founders
        colnames(error_matrix) <- founders
      }
      
      if (verbose >= 3) {
        show_clustering_diagnostics(distances, founders, h_cutoff, groups, n_groups, haplotype_freqs)
      }
      
      estimate_OK <- ifelse(n_groups == length(founders), 1, 0)
      
    } else {
      # LSEI failed with error code
      estimate_OK <- NA
      # Only show diagnostics on actual error
      cat("✗ LSEI failed with error code:", lsei_result$IsError, "\n")
      if (lsei_result$IsError == 3) {
        cat("  Inequalities are contradictory - showing diagnostics:\n")
        # Show founder ranges for this problematic case
        for (i in seq_along(founders)) {
          freq_range <- range(founder_matrix_clean[, i], na.rm = TRUE)
          cat(sprintf("  %s: %.3f - %.3f\n", founders[i], freq_range[1], freq_range[2]))
        }
        # Show matrix condition if available
        tryCatch({
          cond <- kappa(founder_matrix_clean)
          cat("  Matrix condition number:", format(cond, scientific = TRUE), "\n")
        }, error = function(e) {})
      }
    }
  }, error = function(e) {
    # Catastrophic LSEI error
    estimate_OK <- NA
    if (verbose >= 1) {
      cat("✗ LSEI error:", e$message, "\n")
      cat("  Matrix dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
      cat("  Complete SNPs:", sum(complete.cases(founder_matrix_clean)), "\n")
    }
  })
  
  return(list(estimate_OK = estimate_OK, haplotype_freqs = haplotype_freqs, 
              groups = groups, error_matrix = error_matrix))
}

#' Create empty result for failed cases
create_empty_result <- function(chr, pos, sample_name, method, window_size, founders) {
  result <- list(
    chr = chr,
    pos = pos,
    sample = sample_name,
    method = method,
    final_window_size = window_size,
    n_snps = 0,
    estimate_OK = NA
  )
  
  # Add founder frequencies as NA
  for (founder in founders) {
    result[[founder]] <- NA
  }
  
  return(result)
}

#' Create result structure
create_result <- function(chr, pos, sample_name, method, window_size, n_snps, estimate_OK, haplotype_freqs, founders, h_cutoff = NA, groups = NULL, error_matrix = NULL) {
  result <- list(
    chr = chr,
    pos = pos,
    sample = sample_name,
    method = method
  )
  
  # Add h_cutoff for adaptive method (right after method)
  if (method == "adaptive") {
    result[["h_cutoff"]] <- h_cutoff
  }
  
  # Continue with remaining columns
  result[["final_window_size"]] <- window_size
  result[["n_snps"]] <- n_snps
  result[["estimate_OK"]] <- estimate_OK
  
  # Add founder frequencies
  for (i in seq_along(founders)) {
    result[[founders[i]]] <- haplotype_freqs[i]
  }
  
  # Add groups and error matrix if provided
  if (!is.null(groups)) {
    result[["groups"]] <- groups
  }
  if (!is.null(error_matrix)) {
    result[["error_matrix"]] <- error_matrix
  }
  
  return(result)
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
