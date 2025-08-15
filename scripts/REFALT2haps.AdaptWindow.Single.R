#!/usr/bin/env Rscript

# REFALT2haps Adaptive Window - Single Parameter Version
# This script runs haplotype estimation for a single adaptive h_cutoff value
# 
# Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>
# Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 2.5

library(tidyverse)
library(limSolve)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript REFALT2haps.AdaptWindow.Single.R <chr> <parfile> <mydir> <h_cutoff>\n")
  cat("Example: Rscript REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/haplotype_parameters.R process/test 2.5\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
h_cutoff_value <- as.numeric(args[4])

# Source the parameter file
source(parfile)

# Override the h_cutoff from parameter file with the provided value
h_cutoff <- h_cutoff_value

# Define file paths
# REFALT files are in the parent directory, not the results subdirectory
filein <- paste0(dirname(mydir), "/RefAlt.", mychr, ".txt")

cat("=== REFALT2haps Adaptive Window - Single Parameter ===\n")
cat("Chromosome:", mychr, "\n")
cat("H_cutoff:", h_cutoff, "\n")
cat("Input file:", filein, "\n\n")

# Load data
cat("Loading data...\n")
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

# Subset dataset to only high-quality SNPs
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get all non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("Testing h_cutoff:", h_cutoff, "\n")
cat("Non-founder samples:", length(non_founder_samples), "\n")
cat("Samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Define scanning positions (500kb to end-500kb, 10kb steps)
chromosome_length <- max(df$POS)
scan_start <- 500000
scan_end <- chromosome_length - 500000
scan_positions <- seq(scan_start, scan_end, by = 10000)

cat("Chromosome length:", chromosome_length, "bp\n")
cat("Scanning from:", scan_start, "to", scan_end, "bp\n")
cat("Total positions to scan:", length(scan_positions), "\n\n")

# Initialize results table
results_list <- list()

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

# Scan each position
for (pos_idx in seq_along(scan_positions)) {
  test_pos <- scan_positions[pos_idx]
  
  if (pos_idx %% 100 == 0) {
    cat("Processing position", pos_idx, "of", length(scan_positions), "(", test_pos, "bp)\n")
  }
  
  # Process each sample
  for (sample_name in non_founder_samples) {
    # Get sample data
    sample_data <- df3 %>%
      filter(name == sample_name) %>%
      select(POS, freq, N)
    
    # Check if we have enough data for estimation
    if (nrow(sample_data) < 10) {
      # Return NA for insufficient sample data
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = NA,
        n_snps = nrow(sample_data)
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Adaptive window haplotype estimation using fixed algorithm
    # Test multiple window sizes with hierarchical clustering and constraint accumulation
    
    # Define window sizes to try (progressive expansion)
    base_window_size <- 10000  # Start with 10kb
    window_sizes <- c(base_window_size, 
                      base_window_size * 2.5, 
                      base_window_size * 5, 
                      base_window_size * 10, 
                      base_window_size * 20,
                      base_window_size * 50)  # Max 500kb
    
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
        filter(CHROM == mychr &
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
        next  # Skip to next window size
      }
      
      # Convert to matrix for clustering
      founder_matrix <- as.matrix(founder_matrix)
      
      # Cluster founders based on similarity using hierarchical clustering
      founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
      
      # Check if groups meaningfully changed
      groups_meaningfully_changed <- groups_changed(founder_clusters, previous_groups)
      
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
        
        # Get unique clusters
        unique_clusters <- unique(founder_clusters)
        
        # Solve constrained least squares
        tryCatch({
          result <- limSolve::lsei(A = founder_matrix, B = sample_freqs, E = E, F = F, 
                                  G = diag(n_founders), H = matrix(rep(0.0003, n_founders)))
          
          if (result$IsError == 0) {
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
              # If no constraints (e.g., all founders in 1 group), reset to NULL
              accumulated_constraints <- NULL
              accumulated_constraint_values <- NULL
            }
            
            # Store the best result for this position/h_cutoff
            best_result <- result
            best_n_groups <- length(unique(founder_clusters))
            
            # Check if all founders are separated
            if (best_n_groups == length(founders)) {
              break  # Stop for this h_cutoff
            }
          }
        }, error = function(e) {
          # Continue to next window size
        })
      }
      
      # Update previous groups for next iteration
      previous_groups <- founder_clusters
    }
    
    # Use the best result found
    if (!is.null(best_result)) {
      founder_frequencies <- best_result$X
    } else {
      # Return NA if no successful estimation
      founder_frequencies <- rep(NA, length(founders))
    }
    
    # Store results
    result_row <- list(
      chr = mychr,
      pos = test_pos,
      sample = sample_name,
      h_cutoff = h_cutoff,
      n_groups = best_n_groups,
      n_snps = nrow(sample_data)
    )
    
    # Add founder frequencies as named columns
    for (i in seq_along(founders)) {
      result_row[[founders[i]]] <- founder_frequencies[i]
    }
    
    results_list[[length(results_list) + 1]] <- result_row
  }
}

# Convert results to data frame
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  # Save results
  output_file <- paste0(mydir, "/adaptive_window_h", h_cutoff, "_results_", mychr, ".RDS")
  saveRDS(results_df, output_file)
  
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("Total haplotype estimates:", nrow(results_df), "\n")
  cat("Output file:", output_file, "\n")
  cat("File size:", file.size(output_file), "bytes\n")
  
  # Calculate success rate (only count non-NA estimates as successful)
  total_positions <- length(scan_positions) * length(non_founder_samples)
  
  # Count successful estimates (where at least one founder frequency is not NA)
  # Check the first founder column (B1) to determine if estimation was successful
  successful_estimates <- results_df %>%
    filter(!is.na(B1)) %>%
    nrow()
  
  success_rate <- successful_estimates / total_positions * 100
  
  cat("\nSuccess rate:", round(success_rate, 1), "% (", successful_estimates, "of", total_positions, "position/sample combinations)\n")
  cat("Total results (including NA):", nrow(results_df), "\n")
  cat("Positions scanned:", length(scan_positions), "\n")
  cat("Samples processed:", length(non_founder_samples), "\n")
  
  # Show sample summary with success rates
  cat("\nEstimates per sample:\n")
  sample_counts <- results_df %>%
    group_by(sample) %>%
    summarize(
      total_results = n(),
      successful_estimates = sum(!is.na(B1)),
      success_rate = successful_estimates / length(scan_positions) * 100
    ) %>%
    arrange(desc(successful_estimates))
  print(sample_counts)
  
  # Show group summary
  cat("\nAverage number of founder groups per estimate:\n")
  cat("Mean:", mean(results_df$n_groups, na.rm = TRUE), "\n")
  cat("Range:", range(results_df$n_groups, na.rm = TRUE), "\n")
  
} else {
  cat("\n❌ No haplotype estimates obtained!\n")
  cat("This could be due to:\n")
  cat("- H_cutoff too low (all founders grouped together)\n")
  cat("- H_cutoff too high (no founder grouping)\n")
  cat("- Insufficient SNP coverage\n")
  cat("- Data quality issues\n")
}
