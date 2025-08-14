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
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")

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
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Get founder data
    founder_data <- df3 %>%
      filter(name %in% founders) %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    if (ncol(founder_data) < length(founders) + 2) {
      # Return NA for missing founders
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = NA,
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Prepare founder matrix (exclude POS column)
    founder_matrix <- founder_data %>%
      select(-POS) %>%
      as.matrix()
    
    # Remove rows with any NA values
    complete_rows <- complete.cases(founder_matrix)
    if (sum(complete_rows) < 10) {
      # Return NA for insufficient complete rows
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = NA,
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    founder_matrix <- founder_matrix[complete_rows, ]
    sample_freqs <- sample_data$freq[match(founder_data$POS[complete_rows], sample_data$POS)]
    
    # Remove NA sample frequencies
    valid_indices <- !is.na(sample_freqs)
    if (sum(valid_indices) < 10) {
      # Return NA for insufficient valid frequencies
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = NA,
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    founder_matrix <- founder_matrix[valid_indices, ]
    sample_freqs <- sample_freqs[valid_indices]
    
    # Calculate pairwise distances between founders
    founder_distances <- as.matrix(dist(t(founder_matrix)))
    
    # Find founder groups based on h_cutoff
    founder_groups <- list()
    used_founders <- rep(FALSE, ncol(founder_matrix))
    
    for (i in 1:ncol(founder_matrix)) {
      if (used_founders[i]) next
      
      # Find all founders within h_cutoff distance
      group_members <- which(founder_distances[i, ] <= h_cutoff & !used_founders)
      if (length(group_members) > 0) {
        founder_groups[[length(founder_groups) + 1]] <- group_members
        used_founders[group_members] <- TRUE
      }
    }
    
    # If no groups found or only one group, return NA
    if (length(founder_groups) <= 1) {
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = length(founder_groups),
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      next
    }
    
    # Create reduced founder matrix by averaging within groups
    reduced_founder_matrix <- matrix(0, nrow = nrow(founder_matrix), ncol = length(founder_groups))
    group_names <- character(length(founder_groups))
    
    for (g in seq_along(founder_groups)) {
      group_founders <- founder_groups[[g]]
      reduced_founder_matrix[, g] <- rowMeans(founder_matrix[, group_founders, drop = FALSE])
      group_names[g] <- paste0("Group", g, "_", paste(founders[group_founders], collapse = "_"))
    }
    
    # Solve constrained least squares with reduced founder matrix
    tryCatch({
      result <- limSolve::lsei(A = reduced_founder_matrix, B = sample_freqs, 
                              E = matrix(1, nrow = 1, ncol = ncol(reduced_founder_matrix)), 
                              F = 1, G = diag(ncol(reduced_founder_matrix)), H = rep(0, ncol(reduced_founder_matrix)))
      
      if (result$IsError == 0) {
        # Expand group frequencies back to individual founders
        founder_frequencies <- rep(0, length(founders))
        for (g in seq_along(founder_groups)) {
          group_founders <- founder_groups[[g]]
          founder_frequencies[group_founders] <- result$X[g] / length(group_founders)
        }
        
        # Store results
        result_row <- list(
          chr = mychr,
          pos = test_pos,
          sample = sample_name,
          h_cutoff = h_cutoff,
          n_groups = length(founder_groups),
          n_snps = nrow(sample_data),
          founder_frequencies = founder_frequencies
        )
        
        # Add founder frequencies as named columns
        for (i in seq_along(founders)) {
          result_row[[founders[i]]] <- founder_frequencies[i]
        }
        
        results_list[[length(results_list) + 1]] <- result_row
      } else {
        # Return NA for lsei error
        result_row <- list(
          chr = mychr,
          pos = test_pos,
          sample = sample_name,
          h_cutoff = h_cutoff,
          n_groups = length(founder_groups),
          n_snps = nrow(sample_data),
          founder_frequencies = rep(NA, length(founders))
        )
        
        # Add founder frequencies as named columns (all NA)
        for (i in seq_along(founders)) {
          result_row[[founders[i]]] <- NA
        }
        
        results_list[[length(results_list) + 1]] <- result_row
      }
    }, error = function(e) {
      # Return NA for any error
      result_row <- list(
        chr = mychr,
        pos = test_pos,
        sample = sample_name,
        h_cutoff = h_cutoff,
        n_groups = length(founder_groups),
        n_snps = nrow(sample_data),
        founder_frequencies = rep(NA, length(founders))
      )
      
      # Add founder frequencies as named columns (all NA)
      for (i in seq_along(founders)) {
        result_row[[founders[i]]] <- NA
      }
      
      results_list[[length(results_list) + 1]] <- result_row
    })
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
  
  # Calculate success rate
  total_positions <- length(scan_positions) * length(non_founder_samples)
  success_rate <- nrow(results_df) / total_positions * 100
  
  cat("\nSuccess rate:", round(success_rate, 1), "% (", nrow(results_df), "of", total_positions, "position/sample combinations)\n")
  cat("Positions scanned:", length(scan_positions), "\n")
  cat("Samples processed:", length(non_founder_samples), "\n")
  
  # Show sample summary with success rates
  cat("\nEstimates per sample:\n")
  sample_counts <- results_df %>%
    group_by(sample) %>%
    summarize(
      n_estimates = n(),
      success_rate = n() / length(scan_positions) * 100
    ) %>%
    arrange(desc(n_estimates))
  print(sample_counts)
  
  # Show group summary
  cat("\nAverage number of founder groups per estimate:\n")
  cat("Mean:", mean(results_df$n_groups), "\n")
  cat("Range:", range(results_df$n_groups), "\n")
  
} else {
  cat("\n❌ No haplotype estimates obtained!\n")
  cat("This could be due to:\n")
  cat("- H_cutoff too low (all founders grouped together)\n")
  cat("- H_cutoff too high (no founder grouping)\n")
  cat("- Insufficient SNP coverage\n")
  cat("- Data quality issues\n")
}
