#!/usr/bin/env Rscript

# =============================================================================
# REFALT2haps Parameter Assessment Script
# =============================================================================
# This script assesses how window size and h_cutoff affect haplotype estimation
# at a specific genomic midpoint without adaptive constraints
# Usage: Rscript scripts/REFALT2haps.Assessment.R chr parfile mydir start_pos end_pos

library(tidyverse)
library(lazy_dt)
library(limSolve)
library(patchwork)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript REFALT2haps.Assessment.R chr parfile mydir start_pos end_pos\n")
  cat("Example: Rscript REFALT2haps.Assessment.R chr2L helpfiles/haplotype_parameters.R process/test 18000000 20000000\n")
  quit(status = 1)
}

mychr <- args[1]
parfile <- args[2]
mydir <- args[3]
start_pos <- as.numeric(args[4])
end_pos <- as.numeric(args[5])

# Calculate midpoint
midpoint <- (start_pos + end_pos) / 2
window_size <- end_pos - start_pos

# Source the parameter file
source(parfile)

# Define file paths
filein <- paste0(mydir, "/RefAlt.", mychr, ".txt")
rdsfile <- paste0(mydir, "/df3.", mychr, ".RDS")

cat("=== REFALT2haps Parameter Assessment ===\n")
cat("Chromosome:", mychr, "\n")
cat("Region:", start_pos, "-", end_pos, "bp\n")
cat("Midpoint:", midpoint, "bp\n")
cat("Window size:", window_size, "bp\n")
cat("Input file:", filein, "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Estimate haplotype frequencies for a single sample at a specific position
#' @param sample_data SNP data for one sample
#' @param target_pos Target position for estimation
#' @param window_size Window size around target position
#' @param h_cutoff Hierarchical clustering cutoff
#' @return List with haplotype frequencies and metadata
estimate_haplotype_at_position <- function(sample_data, target_pos, window_size, h_cutoff) {
  # Create window around target position
  window_start <- max(0, target_pos - window_size/2)
  window_end <- target_pos + window_size/2
  
  # Filter SNPs in window
  window_snps <- sample_data %>%
    filter(pos >= window_start, pos <= window_end) %>%
    filter(!is.na(freq))
  
  if (nrow(window_snps) == 0) {
    return(list(
      frequencies = rep(NA, length(founders)),
      n_snps = 0,
      window_start = window_start,
      window_end = window_end,
      h_cutoff = h_cutoff
    ))
  }
  
  # Extract founder matrix and sample frequencies
  founder_matrix <- window_snps %>% select(matches(founders))
  sample_freqs <- window_snps$freq
  
  # Convert to matrix for clustering
  founder_matrix <- as.matrix(founder_matrix)
  
  # Cluster founders based on similarity
  founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = h_cutoff)
  
  # Build constraint matrix
  n_founders <- ncol(founder_matrix)
  E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
  F <- 1.0
  
  # Add group constraints for each cluster
  unique_clusters <- unique(founder_clusters)
  for (cluster_id in unique_clusters) {
    cluster_founders <- which(founder_clusters == cluster_id)
    if (length(cluster_founders) > 1) {
      constraint_row <- rep(0, n_founders)
      constraint_row[cluster_founders] <- 1
      E <- rbind(E, constraint_row)
      F <- c(F, 1.0)  # Each group sums to 1
    }
  }
  
  # Solve constrained least squares
  tryCatch({
    result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F)
    
    list(
      frequencies = result$X,
      n_snps = nrow(window_snps),
      window_start = window_start,
      window_end = window_end,
      h_cutoff = h_cutoff,
      error = result$residualNorm
    )
  }, error = function(e) {
    list(
      frequencies = rep(NA, length(founders)),
      n_snps = nrow(window_snps),
      window_start = window_start,
      window_end = window_end,
      h_cutoff = h_cutoff,
      error = NA
    )
  })
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# Load data
cat("Loading data...\n")
if (!file.exists(filein)) {
  stop("Input file not found: ", filein)
}

# Load the processed SNP data
if (file.exists(rdsfile)) {
  df3 <- readRDS(rdsfile)
} else {
  # Fallback: load from RefAlt file if RDS doesn't exist
  cat("RDS file not found, loading from RefAlt file...\n")
  df3 <- read.table(filein, header = TRUE, sep = "\t")
  colnames(df3)[1:2] <- c("chr", "pos")
}

# Define window sizes to test (in bp) - limited to stay within test region
max_window <- min(1000000, window_size / 2)  # Don't exceed half the test region
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000, max_window)
h_cutoffs <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)

cat("Testing", length(window_sizes), "window sizes (max:", max_window/1000, "kb) and", length(h_cutoffs), "h_cutoff values\n")
cat("Test region:", start_pos, "-", end_pos, "bp (midpoint:", midpoint, "bp)\n")
cat("Maximum window size limited to:", max_window/1000, "kb to stay within region bounds\n\n")

# Get sample names (excluding chr and pos columns)
sample_cols <- colnames(df3)[!colnames(df3) %in% c("chr", "pos")]
if (length(sample_cols) == 0) {
  stop("No sample columns found in data")
}

# Test first sample for demonstration
test_sample <- sample_cols[1]
cat("Testing with sample:", test_sample, "\n")

# Create sample data
sample_data <- df3 %>%
  select(chr, pos, all_of(founders)) %>%
  mutate(freq = df3[[test_sample]])

# =============================================================================
# WINDOW SIZE ANALYSIS (Fixed h_cutoff)
# =============================================================================

cat("Analyzing window size effects using default h_cutoff:", h_cutoff, "\n")
window_results <- list()

for (ws in window_sizes) {
  result <- estimate_haplotype_at_position(sample_data, midpoint, ws, h_cutoff)
  window_results[[as.character(ws)]] <- result
}

# Create window size analysis dataframe
window_df <- map_dfr(window_results, function(result) {
  tibble(
    window_size = result$window_size,
    n_snps = result$n_snps,
    error = result$error,
    founder = founders,
    frequency = result$frequencies
  )
}, .id = "window_size_kb") %>%
  mutate(
    window_size_kb = as.numeric(window_size_kb) / 1000,
    window_size = as.numeric(window_size)
  )

# =============================================================================
# H_CUTOFF ANALYSIS (Fixed Window Size)
# =============================================================================

cat("Analyzing h_cutoff effects using fixed 500kb window\n")
hcutoff_results <- list()

for (hc in h_cutoffs) {
  result <- estimate_haplotype_at_position(sample_data, midpoint, 500000, hc)  # Use 500kb window
  hcutoff_results[[as.character(hc)]] <- result
}

# Create h_cutoff analysis dataframe
hcutoff_df <- map_dfr(hcutoff_results, function(result) {
  tibble(
    h_cutoff = result$h_cutoff,
    n_snps = result$n_snps,
    error = result$error,
    founder = founders,
    frequency = result$frequencies
  )
}, .id = "h_cutoff") %>%
  mutate(h_cutoff = as.numeric(h_cutoff))

# =============================================================================
# ADAPTIVE WINDOW ANALYSIS (Different h_cutoffs)
# =============================================================================

cat("Analyzing adaptive window algorithm with different h_cutoffs\n")
adaptive_results <- list()

# Test adaptive algorithm with different h_cutoffs
for (hc in h_cutoffs) {
  # For adaptive analysis, we test how h_cutoff affects clustering
  # We'll use a reasonable window size as a starting point, but the algorithm
  # will adapt based on the clustering results
  test_window_size <- 500000  # 500kb starting window
  
  # Create window around target position
  window_start <- max(0, midpoint - test_window_size/2)
  window_end <- midpoint + test_window_size/2
  
  # Filter SNPs in window
  window_snps <- sample_data %>%
    filter(pos >= window_start, pos <= window_end) %>%
    filter(!is.na(freq))
  
  if (nrow(window_snps) > 0) {
    # Extract founder matrix and sample frequencies
    founder_matrix <- window_snps %>% select(matches(founders))
    sample_freqs <- window_snps$freq
    
    # Convert to matrix for clustering
    founder_matrix <- as.matrix(founder_matrix)
    
    # Use the current h_cutoff for clustering (this is what we're testing)
    founder_clusters <- cutree(hclust(dist(t(founder_matrix))), h = hc)
    
    # Count founder groups and analyze clustering structure
    n_groups <- length(unique(founder_clusters))
    
    # Analyze group sizes
    group_sizes <- table(founder_clusters)
    max_group_size <- max(group_sizes)
    min_group_size <- min(group_sizes)
    
    # Build constraint matrix based on clustering
    n_founders <- ncol(founder_matrix)
    E <- matrix(rep(1, n_founders), nrow = 1)  # Sum to 1 constraint
    F <- 1.0
    
    # Add group constraints for each cluster
    unique_clusters <- unique(founder_clusters)
    for (cluster_id in unique_clusters) {
      cluster_founders <- which(founder_clusters == cluster_id)
      if (length(cluster_founders) > 1) {
        constraint_row <- rep(0, n_founders)
        constraint_row[cluster_founders] <- 1
        E <- rbind(E, constraint_row)
        F <- c(F, 1.0)  # Each group sums to 1
      }
    }
    
    # Solve constrained least squares
    tryCatch({
      result <- lsei(A = founder_matrix, B = sample_freqs, E = E, F = F)
      
      adaptive_results[[as.character(hc)]] <- list(
        frequencies = result$X,
        n_snps = nrow(window_snps),
        n_groups = n_groups,
        max_group_size = max_group_size,
        min_group_size = min_group_size,
        h_cutoff = hc,
        error = result$residualNorm,
        founder_groups = founder_clusters
      )
    }, error = function(e) {
      adaptive_results[[as.character(hc)]] <- list(
        frequencies = rep(NA, length(founders)),
        n_snps = nrow(window_snps),
        n_groups = n_groups,
        max_group_size = max_group_size,
        min_group_size = min_group_size,
        h_cutoff = hc,
        error = NA,
        founder_groups = founder_clusters
      )
    })
  } else {
    adaptive_results[[as.character(hc)]] <- list(
      frequencies = rep(NA, length(founders)),
      n_snps = 0,
      n_groups = 0,
      max_group_size = 0,
      min_group_size = 0,
      h_cutoff = hc,
      error = NA,
      founder_groups = NULL
    )
  }
}

# Create adaptive analysis dataframe
adaptive_df <- map_dfr(adaptive_results, function(result) {
  tibble(
    h_cutoff = result$h_cutoff,
    n_snps = result$n_snps,
    n_groups = result$n_groups,
    max_group_size = result$max_group_size,
    min_group_size = result$min_group_size,
    error = result$error,
    founder = founders,
    frequency = result$frequencies
  )
}, .id = "h_cutoff") %>%
  mutate(h_cutoff = as.numeric(h_cutoff))

# =============================================================================
# CREATE PLOTS
# =============================================================================

cat("Creating assessment plots...\n")

# Colorblind-friendly palette
color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#990099")

# Plot 1: Window size effects on founder frequencies
p1 <- ggplot(window_df, aes(x = window_size_kb, y = frequency, color = founder)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(
    title = paste("Haplotype Frequency Estimates vs Window Size"),
    subtitle = paste("Position:", format(midpoint, scientific = FALSE), "bp on", mychr, "- Fixed h_cutoff =", h_cutoff),
    x = "Window Size (kb)",
    y = "Estimated Frequency",
    color = "Founder"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Plot 2: H_cutoff effects on founder frequencies
p2 <- ggplot(hcutoff_df, aes(x = h_cutoff, y = frequency, color = founder)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(
    title = paste("Haplotype Frequency Estimates vs H_cutoff"),
    subtitle = paste("Position:", format(midpoint, scientific = FALSE), "bp, Fixed 500kb window"),
    x = "H_cutoff",
    y = "Estimated Frequency",
    color = "Founder"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Plot 3: Number of SNPs vs window size
p3 <- ggplot(window_df %>% distinct(window_size_kb, n_snps), 
             aes(x = window_size_kb, y = n_snps)) +
  geom_line(size = 1.5, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  labs(
    title = "Number of SNPs vs Window Size",
    subtitle = paste("Position:", format(midpoint, scientific = FALSE), "bp on", mychr),
    x = "Window Size (kb)",
    y = "Number of SNPs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Plot 4: Adaptive analysis - Founder frequencies vs h_cutoff
p4 <- ggplot(adaptive_df, aes(x = h_cutoff, y = frequency, color = founder)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(
    title = paste("Adaptive Algorithm: Founder Frequencies vs H_cutoff"),
    subtitle = paste("Position:", format(midpoint, scientific = FALSE), "bp - Testing clustering sensitivity"),
    x = "H_cutoff",
    y = "Estimated Frequency",
    color = "Founder"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Plot 5: Clustering behavior vs h_cutoff
p5 <- ggplot(adaptive_df %>% distinct(h_cutoff, n_groups, max_group_size, min_group_size), 
             aes(x = h_cutoff)) +
  geom_line(aes(y = n_groups, color = "Number of Groups"), size = 1.5) +
  geom_point(aes(y = n_groups, color = "Number of Groups"), size = 3) +
  geom_line(aes(y = max_group_size, color = "Max Group Size"), size = 1.5) +
  geom_point(aes(y = max_group_size, color = "Max Group Size"), size = 3) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(
    title = paste("Adaptive Algorithm: Clustering Behavior vs H_cutoff"),
    subtitle = paste("Position:", format(midpoint, scientific = FALSE), "bp - Group structure analysis"),
    x = "H_cutoff",
    y = "Count/Size",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Combine plots in a 2x3 grid (with last row having only 2 plots)
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save plots
output_file <- paste0(mydir, "/comprehensive_assessment_", mychr, "_", format(midpoint, scientific = FALSE), ".png")
ggsave(output_file, combined_plot, width = 16, height = 12, dpi = 300)

# Save data
output_data <- paste0(mydir, "/comprehensive_assessment_", mychr, "_", format(midpoint, scientific = FALSE), ".RDS")
saveRDS(list(
  window_analysis = window_df,
  hcutoff_analysis = hcutoff_df,
  adaptive_analysis = adaptive_df,
  parameters = list(
    chromosome = mychr,
    midpoint = midpoint,
    start_pos = start_pos,
    end_pos = end_pos,
    sample = test_sample
  )
), output_data)

cat("\n=== Comprehensive Assessment Complete ===\n")
cat("Plots saved to:", output_file, "\n")
cat("Data saved to:", output_data, "\n")
cat("Tested window sizes:", paste(window_sizes/1000, "kb"), "\n")
cat("Tested h_cutoffs:", paste(h_cutoffs), "\n")
cat("Sample analyzed:", test_sample, "\n")
cat("\nAnalysis includes:\n")
cat("- Fixed window size effects (using default h_cutoff)\n")
cat("- Fixed h_cutoff effects (using 500kb window)\n")
cat("- Adaptive algorithm effects (using different h_cutoffs)\n")
