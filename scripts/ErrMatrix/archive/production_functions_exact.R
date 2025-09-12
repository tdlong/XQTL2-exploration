# EXACT COPY of production functions - NO MODIFICATIONS
# Only the functions needed for adaptive h4 estimation

library(tidyverse)
library(limSolve)
library(MASS)

# =============================================================================
# EXACT COPY OF estimate_haplotypes_list_format FROM PRODUCTION
# =============================================================================

estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  # CRITICAL: This function's return format is HARDWIRED and CANNOT be changed
  # The downstream pipeline expects exactly this structure:
  # - Groups: integer vector of cluster assignments
  # - Haps: named numeric vector of founder frequencies  
  # - Err: matrix of error/covariance estimates
  # - Names: character vector of founder names
  # Any changes to this format will break the entire pipeline
  
  # No redundant debug output - function is called for every position
  
  method <- match.arg(method)
  
  # No redundant position processing output - already shown by wrapper
    # ADAPTIVE WINDOW METHOD
    if (verbose >= 2) {
      cat(sprintf("=== ADAPTIVE WINDOW: h_cutoff=%g, pos=%s, sample=%s ===\n", 
                  h_cutoff, format(pos, big.mark=","), sample_name))
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
        cat(sprintf("  Window %d: %s", window_idx, 
                    ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp"))))
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
        cat(sprintf(" - %d SNPs, %d groups", nrow(founder_matrix_clean), n_groups))
      }
      
      # Check if clustering improved
      if (window_idx > 1 && n_groups <= previous_n_groups) {
        if (verbose >= 2) cat(" - ✗ No improvement\n")
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
            cat(sprintf(" - ✓ LSEI success, %d groups", n_groups))
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
              
              # Constraint building details removed for cleaner output
            } else {
              # Single founder: lock their exact frequency
              founder_freq <- result$X[cluster_founders]
              
              # Create constraint: this founder = their exact frequency
              constraint_row <- rep(0, n_founders)
              constraint_row[cluster_founders] <- 1
              
              current_constraints <- rbind(current_constraints, constraint_row)
              current_constraint_values <- c(current_constraint_values, founder_freq)
              
              # Constraint building details removed for cleaner output
            }
          }
          
          # Update accumulated constraints for next window
          if (!is.null(current_constraints)) {
            accumulated_constraints <- current_constraints
            accumulated_constraint_values <- current_constraint_values
            if (verbose >= 2) {
              cat(sprintf(" - %d constraints", nrow(current_constraints)))
            }
          } else {
            accumulated_constraints <- NULL
            accumulated_constraint_values <- NULL
          }
          
          # Store the result (will be overwritten as we expand)
          final_result <- result
          final_n_groups <- n_groups
          
          if (verbose >= 2) {
            cat(sprintf("    Stored LSEI result with %d groups\n", n_groups))
            if ("covar" %in% names(result)) {
              cat("    Covariance matrix available in this result (covar field)\n")
            } else {
              cat("    No covariance matrix in this result\n")
            }
          }
          
          # Check if all founders are separated
          if (n_groups == length(founders)) {
            if (verbose >= 2) {
              cat(" - ✅ SUCCESS: All founders distinguished\n")
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

    # Return results
    if (is.null(final_result)) {
      return(NULL)
    }
    
    # Calculate founder frequencies and error matrix
    founder_frequencies <- final_result$X
    names(founder_frequencies) <- founders
    
    # Use covariance matrix if available, otherwise create identity
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      error_matrix <- final_result$covar
    } else {
      # Fallback: create identity matrix
      error_matrix <- diag(length(founders))
    }
    
    # Check for singular matrix
    kappa_val <- kappa(error_matrix)
    if (is.infinite(kappa_val) || kappa_val > 1e12) {
      if (verbose >= 1) {
        cat("Warning: Singular error matrix (kappa =", kappa_val, "), using pseudoinverse\n")
      }
      error_matrix <- ginv(error_matrix)
    }
    
    return(list(
      Groups = groups,
      Haps = founder_frequencies,
      Err = error_matrix,
      Names = founders
    ))
}

# =============================================================================
# EXACT COPY OF process_refalt_data FROM PRODUCTION
# =============================================================================

process_refalt_data <- function(refalt_file, founders) {
  # Load RefAlt data and process into df3 format - EXACT from working code
  cat("Loading RefAlt data from:", refalt_file, "\n")
  
  refalt_data <- read.table(refalt_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Process into df3 format (one row per sample per position) - EXACT from working code
  df3 <- refalt_data %>%
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
    select(-REF, -ALT) %>%
    as_tibble()
  
  # Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
  founder_wide <- df3 %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  quality_filtered_positions <- founder_wide %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    ) %>%
    pull(POS)
  
  # Filter to quality positions and include sample data
  df3 <- df3 %>%
    filter(POS %in% quality_filtered_positions)
  
  cat("✓ Processed", nrow(df3), "rows for", length(unique(df3$name)), "samples\n")
  return(df3)
}
