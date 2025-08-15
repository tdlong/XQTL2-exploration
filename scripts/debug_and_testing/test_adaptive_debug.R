#!/usr/bin/env Rscript

# Debug script to test adaptive window haplotype estimation
# Run for 20 positions with 2 different h_cutoff values, showing intermediate steps

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Load parameter file
source("helpfiles/JUICE/JUICE_haplotype_parameters.R")

# Load data
cat("Loading REFALT data...\n")
df <- read.table("process/JUICE/RefAlt.chr2R.txt", header = TRUE)

# Transform to frequencies
cat("Converting to frequencies...\n")
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

# Subset dataset
df3 <- good_snps %>%
  left_join(df2, multiple = "all")

# Get non-founder samples
all_samples <- unique(df2$name)
non_founder_samples <- all_samples[!all_samples %in% founders]

cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("Non-founder samples:", paste(non_founder_samples, collapse = ", "), "\n\n")

# Test 20 positions in the middle of the chromosome
test_positions <- seq(10000000, 11000000, by = 50000)[1:20]  # 20 positions
test_h_cutoffs <- c(4, 10)  # Test two different cutoffs

cat("Testing 20 positions from 10M to 11M bp\n")
cat("Testing h_cutoff values:", paste(test_h_cutoffs, collapse = ", "), "\n\n")

# Test each h_cutoff
for (h_cutoff in test_h_cutoffs) {
  cat("=", strrep("=", 60), "\n", sep = "")
  cat("TESTING h_cutoff =", h_cutoff, "\n")
  cat("=", strrep("=", 60), "\n\n", sep = "")
  
  # Test first 5 positions in detail
  for (pos_idx in 1:5) {
    test_pos <- test_positions[pos_idx]
    cat("--- Position", pos_idx, ":", test_pos, "bp ---\n")
    
    # Get sample data
    sample_name <- non_founder_samples[1]  # Test first sample
    sample_data <- df3 %>%
      filter(name == sample_name) %>%
      select(POS, freq, N)
    
    cat("Sample data:", nrow(sample_data), "SNPs\n")
    
    # Get founder data
    founder_data <- df3 %>%
      filter(name %in% founders) %>%
      select(POS, name, freq) %>%
      pivot_wider(names_from = name, values_from = freq)
    
    cat("Founder data:", nrow(founder_data), "SNPs\n")
    
    # Prepare founder matrix
    founder_matrix <- founder_data %>%
      select(-POS) %>%
      as.matrix()
    
    # Remove rows with any NA values
    complete_rows <- complete.cases(founder_matrix)
    founder_matrix <- founder_matrix[complete_rows, ]
    
    cat("Complete founder matrix:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
    
    if (nrow(founder_matrix) < 10) {
      cat("⚠️  Insufficient data for estimation\n\n")
      next
    }
    
    # Calculate pairwise distances between founders
    founder_distances <- as.matrix(dist(t(founder_matrix)))
    
    cat("Founder distance matrix:\n")
    print(round(founder_distances, 4))
    
    cat("Distance statistics:\n")
    cat("  Min distance:", min(founder_distances[founder_distances > 0]), "\n")
    cat("  Max distance:", max(founder_distances), "\n")
    cat("  Mean distance:", mean(founder_distances[founder_distances > 0]), "\n")
    
    # Find founder groups based on h_cutoff
    cat("\nFinding founder groups with h_cutoff =", h_cutoff, "...\n")
    founder_groups <- list()
    used_founders <- rep(FALSE, ncol(founder_matrix))
    
    for (i in 1:ncol(founder_matrix)) {
      if (used_founders[i]) next
      
      # Find all founders within h_cutoff distance
      group_members <- which(founder_distances[i, ] <= h_cutoff & !used_founders)
      if (length(group_members) > 0) {
        founder_groups[[length(founder_groups) + 1]] <- group_members
        used_founders[group_members] <- TRUE
        cat("  Group", length(founder_groups), ":", paste(founders[group_members], collapse = ", "), "\n")
      }
    }
    
    cat("Total groups:", length(founder_groups), "\n")
    
    if (length(founder_groups) <= 1) {
      cat("⚠️  Too few groups - estimation will fail\n\n")
      next
    }
    
    # Create reduced founder matrix by averaging within groups
    cat("Creating reduced founder matrix...\n")
    reduced_founder_matrix <- matrix(0, nrow = nrow(founder_matrix), ncol = length(founder_groups))
    group_names <- character(length(founder_groups))
    
    for (g in seq_along(founder_groups)) {
      group_founders <- founder_groups[[g]]
      reduced_founder_matrix[, g] <- rowMeans(founder_matrix[, group_founders, drop = FALSE])
      group_names[g] <- paste0("Group", g, "_", paste(founders[group_founders], collapse = "_"))
      cat("  Group", g, "average:", round(mean(reduced_founder_matrix[, g]), 4), "\n")
    }
    
    # Get sample frequencies for these positions
    sample_freqs <- sample_data$freq[match(founder_data$POS[complete_rows], sample_data$POS)]
    sample_freqs <- sample_freqs[!is.na(sample_freqs)]
    
    if (length(sample_freqs) < 10) {
      cat("⚠️  Insufficient sample frequencies\n\n")
      next
    }
    
    cat("Sample frequencies (first 10):", round(head(sample_freqs, 10), 4), "\n")
    
    # Solve constrained least squares
    cat("Solving constrained least squares...\n")
    tryCatch({
      result <- limSolve::lsei(A = reduced_founder_matrix, B = sample_freqs, 
                              E = matrix(1, nrow = 1, ncol = ncol(reduced_founder_matrix)), 
                              F = 1, G = diag(ncol(reduced_founder_matrix)), H = rep(0, ncol(reduced_founder_matrix)))
      
      if (result$IsError == 0) {
        cat("✓ LSEI successful\n")
        cat("Group frequencies:", round(result$X, 4), "\n")
        
        # Expand group frequencies back to individual founders
        founder_frequencies <- rep(0, length(founders))
        for (g in seq_along(founder_groups)) {
          group_founders <- founder_groups[[g]]
          founder_frequencies[group_founders] <- result$X[g] / length(group_founders)
        }
        
        cat("Founder frequencies:", round(founder_frequencies, 4), "\n")
      } else {
        cat("❌ LSEI failed\n")
      }
    }, error = function(e) {
      cat("❌ LSEI error:", e$message, "\n")
    })
    
    cat("\n")
  }
  
  cat("\n")
}
