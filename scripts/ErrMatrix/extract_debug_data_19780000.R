#!/usr/bin/env Rscript

# Extract debugging data for position 19780000 on chr3R
# This script extracts the raw data that gets passed to the haplotype estimator
# for both adaptive and fixed window methods at the problematic position

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript extract_debug_data_19780000.R <param_file> <output_dir> <position>")
}

param_file <- args[1]
output_dir <- args[2]
test_position <- as.numeric(args[3])

cat("=== EXTRACTING DEBUG DATA FOR POSITION", test_position, "===\n")
cat("Parameter file:", param_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Test position:", test_position, "\n\n")

# Load parameters
source(param_file)
cat("✓ Parameters loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")
cat("H cutoff:", h_cutoff, "\n\n")

# Load observed data
refalt_file <- file.path(output_dir, paste0("RefAlt.chr3R.txt"))
if (!file.exists(refalt_file)) {
  stop("RefAlt file not found:", refalt_file)
}

cat("Loading REFALT data from:", refalt_file, "\n")
refalt_data <- read_tsv(refalt_file, col_names = c("chr", "pos", "name", "ref_count", "alt_count"), 
                       col_types = "cdcdd", show_col_types = FALSE)

cat("✓ REFALT data loaded:", nrow(refalt_data), "rows\n")

# Convert counts to frequencies
observed_data <- refalt_data %>%
  mutate(
    total_count = ref_count + alt_count,
    freq = alt_count / total_count
  ) %>%
  filter(total_count > 0) %>%
  select(chr, pos, name, freq)

cat("✓ Converted to frequencies:", nrow(observed_data), "rows\n")

# Filter for euchromatin
euchromatin_boundaries <- list(
  chr2L = c(1, 23011544),
  chr2R = c(5398184, 24684540),
  chr3L = c(1, 24543557),
  chr3R = c(1, 27905053),
  chr4 = c(1, 1351857),
  chrX = c(1, 22422827)
)

chr <- "chr3R"
if (chr %in% names(euchromatin_boundaries)) {
  boundaries <- euchromatin_boundaries[[chr]]
  observed_data <- observed_data %>%
    filter(pos >= boundaries[1] & pos <= boundaries[2])
  cat("✓ Filtered for euchromatin:", nrow(observed_data), "rows\n")
  cat("Euchromatin region:", boundaries[1], "-", boundaries[2], "bp\n")
}

# Get samples
samples <- unique(observed_data$name)
samples <- samples[!samples %in% founders]
cat("Samples found:", length(samples), "samples\n")
cat("First 5 samples:", paste(head(samples, 5), collapse = ", "), "\n\n")

# Test on first sample
test_sample <- samples[1]
cat("Testing on sample:", test_sample, "\n\n")

# Extract data for different window sizes (adaptive method)
window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
debug_data <- list()

for (window_size in window_sizes) {
  cat("=== EXTRACTING DATA FOR WINDOW SIZE:", window_size, "bp ===\n")
  
  window_start <- test_position - window_size/2
  window_end <- test_position + window_size/2
  
  window_data <- observed_data %>%
    filter(pos >= window_start & pos <= window_end & name %in% c(founders, test_sample))
  
  cat("SNPs in window:", nrow(window_data), "\n")
  cat("Window range:", window_start, "-", window_end, "bp\n")
  
  if (nrow(window_data) == 0) {
    cat("❌ No SNPs found in window\n\n")
    next
  }
  
  # Convert to wide format
  wide_data <- window_data %>%
    select(pos, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  cat("Wide data dimensions:", nrow(wide_data), "x", ncol(wide_data), "\n")
  
  if (!all(c(founders, test_sample) %in% names(wide_data))) {
    cat("❌ Missing founders or sample in wide data\n\n")
    next
  }
  
  # Get founder matrix and sample frequencies
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  sample_freqs <- wide_data %>%
    pull(!!test_sample)
  
  # Remove rows with missing data
  complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  cat("Clean data dimensions:", nrow(founder_matrix_clean), "x", ncol(founder_matrix_clean), "\n")
  cat("Complete rows:", sum(complete_rows), "/", nrow(founder_matrix), "\n")
  
  if (nrow(founder_matrix_clean) < 10) {
    cat("❌ Insufficient SNPs after cleaning:", nrow(founder_matrix_clean), "\n\n")
    next
  }
  
  # Test hierarchical clustering
  founder_dist <- dist(t(founder_matrix_clean))
  hclust_result <- hclust(founder_dist, method = "complete")
  groups <- cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))
  
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  cat("Number of groups:", n_groups, "\n")
  
  # Test LSEI
  n_founders <- ncol(founder_matrix_clean)
  E <- matrix(rep(1, n_founders), nrow = 1)
  F <- 1.0
  
  result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                        E = E, F = F, 
                        G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                        fulloutput = TRUE)
  
  cat("LSEI IsError:", result$IsError, "\n")
  
  if (result$IsError == 0) {
    cat("✓ LSEI successful\n")
    hap_freqs <- result$X
    names(hap_freqs) <- founders
    cat("Haplotype frequencies sum:", round(sum(hap_freqs), 4), "\n")
    
    # Store debug data
    debug_data[[as.character(window_size)]] <- list(
      window_size = window_size,
      window_start = window_start,
      window_end = window_end,
      n_snps = nrow(founder_matrix_clean),
      founder_matrix = founder_matrix_clean,
      sample_freqs = sample_freqs_clean,
      groups = groups,
      n_groups = n_groups,
      hap_freqs = hap_freqs,
      error_matrix = result$cov,
      lsei_error = result$IsError,
      wide_data = wide_data
    )
  } else {
    cat("❌ LSEI failed with error code:", result$IsError, "\n")
  }
  
  cat("\n")
}

# Also extract data for fixed window method (100kb)
cat("=== EXTRACTING FIXED WINDOW DATA (100kb) ===\n")
fixed_window_size <- 100000
window_start <- test_position - fixed_window_size/2
window_end <- test_position + fixed_window_size/2

window_data <- observed_data %>%
  filter(pos >= window_start & pos <= window_end & name %in% c(founders, test_sample))

cat("SNPs in fixed window:", nrow(window_data), "\n")
cat("Fixed window range:", window_start, "-", window_end, "bp\n")

if (nrow(window_data) > 0) {
  wide_data <- window_data %>%
    select(pos, name, freq) %>%
    pivot_wider(names_from = name, values_from = freq)
  
  founder_matrix <- wide_data %>%
    select(all_of(founders)) %>%
    as.matrix()
  sample_freqs <- wide_data %>%
    pull(!!test_sample)
  
  complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  if (nrow(founder_matrix_clean) >= 10) {
    founder_dist <- dist(t(founder_matrix_clean))
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(rep(1, n_founders), nrow = 1)
    F <- 1.0
    
    result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                          E = E, F = F, 
                          G = diag(n_founders), H = matrix(rep(0.0003, n_founders)),
                          fulloutput = TRUE)
    
    if (result$IsError == 0) {
      hap_freqs <- result$X
      names(hap_freqs) <- founders
      
      debug_data[["fixed_100kb"]] <- list(
        window_size = fixed_window_size,
        window_start = window_start,
        window_end = window_end,
        n_snps = nrow(founder_matrix_clean),
        founder_matrix = founder_matrix_clean,
        sample_freqs = sample_freqs_clean,
        groups = groups,
        n_groups = n_groups,
        hap_freqs = hap_freqs,
        error_matrix = result$cov,
        lsei_error = result$IsError,
        wide_data = wide_data
      )
      
      cat("✓ Fixed window data extracted successfully\n")
    } else {
      cat("❌ Fixed window LSEI failed\n")
    }
  }
}

# Save debug data
debug_file <- paste0("debug_data_position_", test_position, ".rds")
saveRDS(debug_data, debug_file)
cat("\n✓ Debug data saved to:", debug_file, "\n")

# Create summary
summary_data <- data.frame(
  method = names(debug_data),
  window_size = sapply(debug_data, function(x) x$window_size),
  n_snps = sapply(debug_data, function(x) x$n_snps),
  n_groups = sapply(debug_data, function(x) x$n_groups),
  lsei_error = sapply(debug_data, function(x) x$lsei_error),
  hap_sum = sapply(debug_data, function(x) round(sum(x$hap_freqs), 4)),
  error_sum = sapply(debug_data, function(x) {
    if (!is.null(x$error_matrix)) round(sum(diag(x$error_matrix)), 6) else NA
  })
)

summary_file <- paste0("debug_summary_position_", test_position, ".csv")
write_csv(summary_data, summary_file)
cat("✓ Summary saved to:", summary_file, "\n")

cat("\n=== EXTRACTION COMPLETE ===\n")
cat("Files created:\n")
cat("-", debug_file, "(RDS file with all debug data)\n")
cat("-", summary_file, "(CSV summary)\n")
cat("\nYou can now download these files for local debugging.\n")
