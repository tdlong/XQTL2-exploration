#!/usr/bin/env Rscript

# Debug script to investigate adaptive window h_cutoff issue
# This script examines founder distances and grouping behavior

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

# Test a few positions to examine founder distances
test_positions <- c(10000000, 15000000, 20000000)  # Test positions

for (test_pos in test_positions) {
  cat("=== Testing position:", test_pos, "===\n")
  
  # Get sample data
  sample_name <- non_founder_samples[1]  # Test first sample
  sample_data <- df3 %>%
    filter(name == sample_name) %>%
    select(POS, freq, N)
  
  # Get founder data
  founder_data <- df3 %>%
    filter(name %in% founders) %>%
    select(POS, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  # Prepare founder matrix
  founder_matrix <- founder_data %>%
    select(-POS) %>%
    as.matrix()
  
  # Remove rows with any NA values
  complete_rows <- complete.cases(founder_matrix)
  founder_matrix <- founder_matrix[complete_rows, ]
  
  cat("Founder matrix dimensions:", nrow(founder_matrix), "x", ncol(founder_matrix), "\n")
  
  # Calculate pairwise distances between founders
  founder_distances <- as.matrix(dist(t(founder_matrix)))
  
  cat("Founder distance matrix:\n")
  print(round(founder_distances, 4))
  
  cat("Distance statistics:\n")
  cat("Min distance:", min(founder_distances[founder_distances > 0]), "\n")
  cat("Max distance:", max(founder_distances), "\n")
  cat("Mean distance:", mean(founder_distances[founder_distances > 0]), "\n")
  cat("Median distance:", median(founder_distances[founder_distances > 0]), "\n")
  
  # Test different h_cutoff values
  h_cutoffs <- c(4, 6, 8, 10)
  
  for (h_cutoff in h_cutoffs) {
    cat("\n--- Testing h_cutoff =", h_cutoff, "---\n")
    
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
    
    cat("Number of groups:", length(founder_groups), "\n")
    
    for (g in seq_along(founder_groups)) {
      group_founders <- founder_groups[[g]]
      cat("  Group", g, ":", paste(founders[group_founders], collapse = ", "), "\n")
    }
    
    if (length(founder_groups) <= 1) {
      cat("  ⚠️  Too few groups - estimation will fail\n")
    }
  }
  
  cat("\n" + strrep("=", 50) + "\n\n")
}
