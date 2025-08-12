#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Adaptive Window Testing Script - Chromosome Scanner
# =============================================================================
# Script to test adaptive window algorithm across entire chromosome
# Usage: Rscript scripts/REFALT2haps.AdaptWindow.R chr parfile mydir [verbose]

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.R chr parfile mydir [verbose]\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.R chr2L helpfiles/haplotype_parameters.R process/test\n")
  cat("Add 'verbose' as 4th argument for detailed output\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
verbose <- ifelse(length(args) >= 4 && args[4] == "verbose", TRUE, FALSE)

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Adaptive Window Test - Chromosome Scanner ===\n")
cat("Chromosome:", mychr, "\n")
cat("Input file:", filein, "\n")
cat("Verbose mode:", verbose, "\n\n")

# Load data
if (verbose) cat("Loading data...\n")
df <- read.table(filein, header = TRUE)

# Transform REF/ALT counts to frequencies
if (verbose) cat("Converting counts to frequencies...\n")
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
if (verbose) cat("Filtering for high-quality SNPs...\n")
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

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Define window sizes to try (progressive expansion, max 500kb as requested)
base_window_size <- 10000  # Start with 10kb
window_sizes <- c(base_window_size, 
                  base_window_size * 2.5, 
                  base_window_size * 5, 
                  base_window_size * 10, 
                  base_window_size * 20,
                  base_window_size * 50)  # Max 500kb

# Define h_cutoff values to test (removed 40 as requested)
h_cutoffs <- c(2, 4, 6, 8, 10, 20)

# Get all non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

if (verbose) {
  cat("Testing h_cutoffs:", paste(h_cutoffs, collapse = ", "), "\n")
  cat("Window sizes:", paste(window_sizes/1000, "kb", collapse = " → "), "\n")
  cat("Non-founder samples:", length(non_founder_samples), "\n")
  cat("Samples:", paste(non_founder_samples, collapse = ", "), "\n\n")
}

# Define scanning positions (500kb to end-500kb, 10kb steps)
chromosome_length <- max(df$POS)
scan_start <- 500000
scan_end <- chromosome_length - 500000
scan_positions <- seq(scan_start, scan_end, by = 10000)

if (verbose) {
  cat("Chromosome length:", chromosome_length, "bp\n")
  cat("Scanning from:", scan_start, "to", scan_end, "bp\n")
  cat("Total positions to scan:", length(scan_positions), "\n\n")
}

# Initialize results table
results_list <- list()

# Scan each position
for (pos_idx in seq_along(scan_positions)) {
  test_pos <- scan_positions[pos_idx]
  
  if (verbose) {
    cat("Position", pos_idx, "/", length(scan_positions), ":", test_pos, "bp\n")
  }
  
  # Test each h_cutoff for this position
  for (hc_idx in seq_along(h_cutoffs)) {
    hc <- h_cutoffs[hc_idx]
    
    if (verbose) {
      cat("  Testing h_cutoff:", hc, "\n")
    }
    
    # Initialize constraints for this h_cutoff
    accumulated_constraints <- NULL
    accumulated_constraint_values <- NULL
    
    # Run adaptive window algorithm for this h_cutoff
    previous_n_groups <- 0  # Track clustering progress
    
    for (window_idx in seq_along(window_sizes)) {
      current_window_size <- window_sizes[window_idx]
      
      if (verbose) {
        cat("    --- Window size:", current_window_size/1000, "kb ---\n")
      }
      
      # Create window around test position
      window_start <- max(0, test_pos - current_window_size/2)
      window_end <- test_pos + current_window_size/2
      
      if (verbose) {
        cat("    Window:", window_start, "-", window_end, "bp (centered at", test_pos, "bp)\n")
      }
      
      # Filter SNPs in this expanding window
      window_snps_current <- df3 %>%
        filter(CHROM == mychr &
               POS > window_start &
               POS < window_end &
               (name %in% founders | name %in% non_founder_samples)) %>%
        select(-c(CHROM, N)) %>%
        pivot_wider(names_from = name, values_from = freq)
      
      # Get non-founder sample columns (exclude POS and founder columns)
      founder_cols <- names(window_snps_current)[names(window_snps_current) %in% founders]
      non_founder_cols <- names(window_snps_current)[!names(window_snps_current) %in% c("POS", founder_cols)]
      
      if (length(non_founder_cols) > 0) {
        window_snps_current <- window_snps_current %>%
          pivot_longer(all_of(non_founder_cols), 
                      names_to = "sample", values_to = "freq") %>%
          select(-POS)
      } else {
        # No non-founder samples found, create empty data frame
        window_snps_current <- data.frame(sample = character(), freq = numeric())
      }
      
      # Count SNPs in this window
      n_snps <- df3 %>%
        filter(CHROM == mychr &
               POS > window_start &
               POS < window_end &
               name %in% founders) %>%
        distinct(POS) %>%
        nrow()
      
      if (verbose) {
        cat("    SNPs in window:", n_snps, "\n")
      }
      
      if (nrow(window_snps_current) == 0) {
        if (verbose) cat("    No SNPs found in window\n")
        next
      }
      
      # Test each non-founder sample
      for (sample_name in non_founder_samples) {
        sample_data_current <- window_snps_current %>% filter(sample == sample_name)
        
        if (nrow(sample_data_current) > 0) {
          # Extract founder matrix and sample frequencies
          founder_matrix <- sample_data_current %>% select(matches(founders))
          sample_freqs <- sample_data_current$freq
          
          # Filter for non-NA values
          valid_positions <- !is.na(sample_freqs)
          sample_freqs <- sample_freqs[valid_positions]
          founder_matrix <- founder_matrix[valid_positions, ]
          
          if (nrow(founder_matrix) > 0) {
            # Convert to matrix for clustering
            founder_matrix <- as.matrix(founder_matrix)
            
            # Cluster founders based on similarity
            founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = hc)
            
            # Count groups
            n_groups <- length(unique(founder_clusters))
            
            if (verbose) {
              cat("    Founder groups:", n_groups, "\n")
            }
            
            # Check if clustering improved (more groups = better separation)
            if (window_idx > 1 && n_groups <= previous_n_groups) {
              if (verbose) {
                cat("    Clustering not improved (", n_groups, " groups vs ", previous_n_groups, "), skipping to next window\n")
              }
              next
            }
            
            # Show group composition
            if (verbose) {
              for (group_id in unique(founder_clusters)) {
                group_founders <- names(founder_clusters[founder_clusters == group_id])
                cat("      Group", group_id, ":", paste(group_founders, collapse = ", "), "\n")
              }
            }
            
            # Update previous_n_groups for next iteration
            previous_n_groups <- n_groups
            
            # Build constraint matrix with accumulated constraints from smaller windows
            n_founders <- ncol(founder_matrix)
            E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
            F <- 1.0
            
            # Add accumulated constraints from previous (smaller) windows
            if (!is.null(accumulated_constraints)) {
              E <- rbind(E, accumulated_constraints)
              F <- c(F, accumulated_constraint_values)
              if (verbose) {
                cat("    Added", nrow(accumulated_constraints), "accumulated constraints\n")
              }
            }
            
            # Define unique_clusters outside the if/else block
            unique_clusters <- unique(founder_clusters)
            
            # For the first window, we don't add group constraints yet
            if (window_idx == 1) {
              if (verbose) cat("    First window: No group constraints yet, just sum to 1\n")
            } else {
              # Add current window group constraints (only for groups with multiple founders)
              # BUT: We can't add constraints without values, so we'll do this after lsei
              multi_founder_groups <- 0
              for (cluster_id in unique_clusters) {
                cluster_founders <- which(founder_clusters == cluster_id)
                if (length(cluster_founders) > 1) {
                  multi_founder_groups <- multi_founder_groups + 1
                }
              }
              if (multi_founder_groups > 0 && verbose) {
                cat("    Will add", multi_founder_groups, "group constraints after lsei\n")
              }
            }
            
            # Solve constrained least squares
            if (verbose) {
              cat("    Constraint matrix E:\n")
              print(E)
              cat("    Constraint values F:\n")
              print(F)
            }
            
            # Sanity checks for bad estimation space
            n_snps <- nrow(founder_matrix)
            n_founders <- ncol(founder_matrix)
            
            # Check 1: Too few SNPs relative to founders (rule of thumb: need at least 3x)
            if (n_snps < n_founders * 3) {
              if (verbose) {
                cat("    ❌ Bad estimation space: ", n_snps, " SNPs for ", n_founders, " founders (need at least ", n_founders * 3, ")\n")
              }
              next  # Skip to next window size
            }
            
            # Check 2: Matrix condition number (numerical stability)
            if (n_snps >= n_founders) {
              condition_num <- kappa(founder_matrix)
              if (condition_num > 1e10) {
                if (verbose) {
                  cat("    ❌ Bad estimation space: Matrix condition number too high (", format(condition_num, scientific = TRUE), ")\n")
                }
                next  # Skip to next window size
              }
            }
            
            # Check 3: Effective rank (how many founders are actually distinguishable)
            if (n_snps >= n_founders) {
              svd_result <- svd(founder_matrix)
              effective_rank <- sum(svd_result$d > 1e-6)
              if (effective_rank < n_founders * 0.7) {
                if (verbose) {
                  cat("    ❌ Bad estimation space: Effective rank too low (", effective_rank, " for ", n_founders, " founders)\n")
                }
                next  # Skip to next window size
              }
            }
            
            if (verbose) {
              cat("    ✓ Estimation space looks good (", n_snps, " SNPs, condition = ", format(kappa(founder_matrix), scientific = TRUE), ")\n")
            }
            
            tryCatch({
              result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                            G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
              
              if (verbose) cat("    ✓ lsei succeeded!\n")
              
              # Show frequency estimates for this window
              if (verbose) {
                cat("    Frequency estimates:\n")
                for (i in seq_along(founders)) {
                  cat("      ", founders[i], ":", sprintf("%.4f", result$X[i]), "\n")
                }
              }
              
              # Accumulate constraints for next (larger) window
              current_constraints <- NULL
              current_constraint_values <- NULL
              
              if (verbose) cat("    Building constraints for next window:\n")
              for (cluster_id in unique_clusters) {
                cluster_founders <- which(founder_clusters == cluster_id)
                if (length(cluster_founders) > 1) {
                  # Create constraint row for this group
                  constraint_row <- rep(0, n_founders)
                  constraint_row[cluster_founders] <- 1
                  
                  # Calculate the actual group frequency from lsei result
                  group_freq <- sum(result$X[cluster_founders])
                  
                  if (verbose) {
                    cat("      Group", cluster_id, "(", paste(names(founder_clusters[founder_clusters == cluster_id]), collapse = ", "), "): ", 
                        sprintf("%.4f", group_freq), " (constraint: sum = ", sprintf("%.4f", group_freq), ")\n")
                  }
                  
                  current_constraints <- rbind(current_constraints, constraint_row)
                  current_constraint_values <- c(current_constraint_values, group_freq)
                } else {
                  # Single founder: lock their exact frequency
                  founder_name <- names(founder_clusters[founder_clusters == cluster_id])
                  founder_freq <- result$X[cluster_founders]
                  
                  if (verbose) {
                    cat("      Group", cluster_id, "(", founder_name, "): ", 
                        sprintf("%.4f", founder_freq), " (locked: ", founder_name, " = ", sprintf("%.4f", founder_freq), ")\n")
                  }
                  
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
                if (verbose) {
                  cat("    ✓ Accumulated", nrow(current_constraints), "group constraints for next window\n")
                }
              } else {
                # If no constraints (e.g., all founders in 1 group), reset to NULL
                accumulated_constraints <- NULL
                accumulated_constraint_values <- NULL
                if (verbose) {
                  cat("    ✓ No meaningful constraints to accumulate, resetting for next window\n")
                }
              }
              
              # Store results for this position/h_cutoff/window combination
              for (founder_idx in seq_along(founders)) {
                founder_name <- founders[founder_idx]
                freq_estimate <- result$X[founder_idx]
                
                results_list[[length(results_list) + 1]] <- list(
                  chr = mychr,
                  pos = test_pos,
                  sample = sample_name,
                  h_cutoff = hc,
                  founder = founder_name,
                  freq = freq_estimate
                )
              }
              
              # Check if all founders are separated
              if (n_groups == length(founders)) {
                if (verbose) {
                  cat("    ✓ All founders separated! Stopping for this h_cutoff.\n")
                }
                break
              }
              
            }, error = function(e) {
              if (verbose) {
                cat("    ✗ Estimation failed:", e$message, "\n")
              }
            })
          }
        }
      }
    }
    
    if (verbose) {
      cat("  === Completed h_cutoff", hc, "===\n")
    }
  }
  
  # Progress indicator
  if (pos_idx %% 100 == 0) {
    cat("Progress:", pos_idx, "/", length(scan_positions), "positions completed\n")
  }
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/adaptive_window_results_", mychr, ".RDS")
  saveRDS(results_df, output_file)
  
  cat("\n=== SCANNING COMPLETE ===\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total estimates:", nrow(results_df), "\n")
  cat("Output format: chr | pos | sample | h_cutoff | founder | freq\n")
  
  # Show summary
  cat("\nSummary by h_cutoff:\n")
  summary_stats <- results_df %>%
    group_by(h_cutoff) %>%
    summarize(
      n_estimates = n(),
      mean_freq = mean(freq, na.rm = TRUE),
      sd_freq = sd(freq, na.rm = TRUE)
    ) %>%
    arrange(h_cutoff)
  
  print(summary_stats)
  
} else {
  cat("\nNo successful estimates obtained!\n")
}
