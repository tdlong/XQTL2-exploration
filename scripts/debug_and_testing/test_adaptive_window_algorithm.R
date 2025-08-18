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
test_h_cutoff <- 10  # Larger h_cutoff to make it more challenging
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

# Run the complete adaptive window algorithm (loop through window_sizes)
cat("Window sizes to test:", paste(window_sizes, collapse = " "), "\n")
cat("=== STARTING ADAPTIVE WINDOW ALGORITHM ===\n")

# Track SNP counts across all window sizes
snps_tracking <- data.frame(
  window_size = window_sizes,
  total_snps = 0,
  after_quality_filter = 0,
  after_founder_filter = 0,
  after_sample_filter = 0,
  final_valid = 0
)

for (window_idx in seq_along(window_sizes)) {
  window_size <- window_sizes[window_idx]
  cat("--- Window", window_idx, ":", window_size, "---\n")
  
  # Calculate window boundaries
  window_start <- test_pos - window_size/2
  window_end <- test_pos + window_size/2
  cat("Window range:", window_start, "to", window_end, "bp\n")
  
  # Get SNPs in this expanding window (FIXED: removed pivot_wider here)
  window_snps <- df3 %>%
    filter(CHROM == "chr2R" &
           POS > window_start &
           POS < window_end &
           (name %in% founders | name == test_sample))
  
  cat("SNPs in window:", nrow(window_snps), "\n")
  
  # Track total SNPs
  snps_tracking$total_snps[window_idx] <- nrow(window_snps)
  
  # Convert to wide format first (rows = positions, columns = sample + founders)
  wide_data <- window_snps %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  cat("Wide format data:", nrow(wide_data), "positions\n")
  
  # Show first 25 SNPs as a table (raw data inspection)
  if (window_idx == 1) {  # Only show for first window to avoid spam
    cat("\n=== FIRST 25 SNPS IN WINDOW (RAW DATA) ===\n")
    
    # Show first 25 positions
    display_data <- wide_data %>%
      arrange(POS) %>%
      head(25) %>%
      mutate(
        POS = as.integer(POS),
        across(-POS, ~round(.x * 100, 0))  # Convert to percentages, no decimals
      )
    
    # Print compact table
    cat("POS     Sample", paste(sprintf("%6s", founders), collapse = ""), "\n")
    cat("----    ------", paste(rep("------", length(founders)), collapse = ""), "\n")
    
    for (i in 1:nrow(display_data)) {
      row <- display_data[i, ]
      sample_col <- test_sample
      cat(sprintf("%6d  %6d", row$POS, row[[sample_col]]))
      for (founder in founders) {
        freq_val <- row[[founder]]
        if (is.na(freq_val)) {
          cat("     --")
        } else {
          cat(sprintf(" %5d", freq_val))
        }
      }
      cat("\n")
    }
    cat("\n")
  }
  
  # Quality filter: keep rows where ALL founders are fixed (freq < 0.03 OR freq > 0.97)
  # Skip positions where any founder has intermediate frequency
  quality_filtered <- wide_data %>%
    filter(
      if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
    )
  
  cat("After quality filter (founders fixed):", nrow(quality_filtered), "positions\n")
  snps_tracking$after_quality_filter[window_idx] <- nrow(quality_filtered)
  
  # Get founder data for this window
  founder_data <- quality_filtered
  
  cat("After founder filter (same as quality):", nrow(founder_data), "positions\n")
  snps_tracking$after_founder_filter[window_idx] <- nrow(founder_data)
  
  # Check if we have enough founder data
  if (ncol(founder_data) < length(founders) + 1) {
    cat("✗ Missing founders, skipping\n")
    next
  }
  
  # Get founder matrix and sample frequencies from wide format
  founder_matrix <- founder_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  
  # Get sample frequencies
  sample_freqs <- founder_data %>%
    pull(!!test_sample)
  
  cat("After sample filter (sample data available):", length(sample_freqs), "positions\n")
  snps_tracking$after_sample_filter[window_idx] <- length(sample_freqs)
  
  # Filter for non-NA values
  valid_positions <- !is.na(sample_freqs)
  sample_freqs <- sample_freqs[valid_positions]
  founder_matrix <- founder_matrix[valid_positions, ]
  
  cat("Final valid positions (no NA):", sum(valid_positions), "out of", length(valid_positions), "\n")
  cat("Sample NA count:", sum(is.na(sample_freqs)), "\n")
  cat("Founder matrix NA count:", sum(is.na(founder_matrix)), "\n")
  snps_tracking$final_valid[window_idx] <- sum(valid_positions)
  
  # Convert to matrix for clustering
  founder_matrix_clean <- founder_matrix[complete.cases(founder_matrix), ]
  
  if (nrow(founder_matrix_clean) < 10) {
    cat("✗ Insufficient data for estimation, skipping\n")
    next
  }
  
  # Debug: Show founder frequency ranges
  cat("Founder frequency ranges:\n")
  for (i in 1:ncol(founder_matrix_clean)) {
    founder_name <- colnames(founder_matrix_clean)[i]
    freq_range <- range(founder_matrix_clean[, i])
    cat(sprintf("  %s: %.4f to %.4f\n", founder_name, freq_range[1], freq_range[2]))
  }
  
  # Hierarchical clustering
  distances <- dist(t(founder_matrix_clean), method = "euclidean")
  hclust_result <- hclust(distances, method = "ward.D2")
  
  # Debug: Show clustering distances
  cat("Clustering distances:\n")
  print(distances)
  
  # Cut tree at h_cutoff
  groups <- cutree(hclust_result, h = test_h_cutoff)
  n_groups <- length(unique(groups))
  
  cat("Founder groups:", n_groups, "\n")
  cat("Group assignments:", paste(groups, collapse=", "), "\n")
  
  # Check if groups meaningfully changed
  groups_meaningfully_changed <- groups_changed(groups, previous_groups)
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
        
        for (cluster_id in unique(groups)) {
          cluster_founders <- which(groups == cluster_id)
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
        best_n_groups <- length(unique(groups))
        
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
  
  previous_groups <- groups
}

# Final results
cat("\n=== ADAPTIVE WINDOW ALGORITHM COMPLETE ===\n")
cat("✓ Algorithm completed successfully\n")
cat("Final number of groups:", best_n_groups, "\n")
if (best_n_groups > 0) {
  cat("Final estimated frequencies:", paste(sprintf("%.4f", best_frequencies), collapse = ", "), "\n")
}

# Print SNP tracking summary
cat("\n=== SNP COUNT PROGRESSION ACROSS WINDOW SIZES ===\n")
print(snps_tracking)

cat("\n=== FILTERING ANALYSIS ===\n")
for (i in 1:nrow(snps_tracking)) {
  if (snps_tracking$total_snps[i] > 0) {
    quality_loss <- snps_tracking$total_snps[i] - snps_tracking$after_quality_filter[i]
    founder_loss <- snps_tracking$after_quality_filter[i] - snps_tracking$after_founder_filter[i]
    sample_loss <- snps_tracking$after_founder_filter[i] - snps_tracking$after_sample_filter[i]
    na_loss <- snps_tracking$after_sample_filter[i] - snps_tracking$final_valid[i]
    
    cat(sprintf("Window %d (%d kb):\n", i, snps_tracking$window_size[i]/1000))
    cat(sprintf("  Quality filter loss: %d SNPs (%.1f%%)\n", quality_loss, 100*quality_loss/snps_tracking$total_snps[i]))
    cat(sprintf("  Founder filter loss: %d SNPs (%.1f%%)\n", founder_loss, 100*founder_loss/snps_tracking$total_snps[i]))
    cat(sprintf("  Sample filter loss: %d SNPs (%.1f%%)\n", sample_loss, 100*sample_loss/snps_tracking$total_snps[i]))
    cat(sprintf("  NA filter loss: %d SNPs (%.1f%%)\n", na_loss, 100*na_loss/snps_tracking$total_snps[i]))
    cat(sprintf("  Final success rate: %.1f%% (%d/%d)\n", 
                100*snps_tracking$final_valid[i]/snps_tracking$total_snps[i],
                snps_tracking$final_valid[i], snps_tracking$total_snps[i]))
  }
}

cat("\n=== FULL PIPELINE TEST COMPLETE ===\n")
