#!/usr/bin/env Rscript

# Clean SNP Imputation - Vectorized Approach
# This follows proper tidyverse principles and is transparent/validatable

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript snp_imputation_clean.R <chr> <param_file> <output_dir> <estimator> <sample_index> [test_n_snps]")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]
sample_index <- as.integer(args[5])
test_n_snps <- if (length(args) >= 6) as.integer(args[6]) else NULL

# Debug mode if test_n_snps is specified
debug_mode <- !is.null(test_n_snps)

cat("=== Clean SNP Imputation ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample index:", sample_index, "\n")
if (debug_mode) {
  cat("DEBUG MODE: Testing with", test_n_snps, "random SNPs\n")
}
cat("\n")

# Load parameters
source(param_file)

# Define founder order consistently
founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Euchromatin boundaries
euchromatin_boundaries <- list(chr2R = c(5398184, 24684540))
euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]

# Load haplotype results
haplotype_file <- if (str_detect(estimator, "^fixed_")) {
  window_size_kb <- str_extract(estimator, "\\d+")
  file.path(output_dir, paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
} else if (str_detect(estimator, "^adaptive_h")) {
  h_cutoff <- str_extract(estimator, "\\d+")
  file.path(output_dir, paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
} else {
  stop("Invalid estimator format")
}

haplotype_results <- read_rds(haplotype_file)
cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Load REFALT data
refalt_file <- file.path(dirname(output_dir), paste0("RefAlt.", chr, ".txt"))
refalt_data <- read_tsv(refalt_file, show_col_types = FALSE)
cat("✓ REFALT data loaded:", nrow(refalt_data), "rows\n")

# Step 1: Process REFALT data to get clean SNP data
cat("Processing REFALT data...\n")
snp_data <- refalt_data %>%
  # Convert to long format
  pivot_longer(
    cols = -c(CHROM, POS),
    names_to = "sample",
    values_to = "count"
  ) %>%
  # Extract REF/ALT and sample name
  mutate(
    ref_alt = str_sub(sample, 1, 3),
    sample_name = str_sub(sample, 5)
  ) %>%
  select(-sample) %>%
  # Convert to wide format with REF/ALT columns
  pivot_wider(
    names_from = ref_alt,
    values_from = count
  ) %>%
  # Calculate frequencies
  mutate(
    freq = REF / (REF + ALT),
    total_count = REF + ALT
  ) %>%
  select(-c(REF, ALT)) %>%
  # Filter for high-quality SNPs
  group_by(CHROM, POS) %>%
  filter(
    # No missing data in founders
    all(total_count > 0 | !(sample_name %in% founders)),
    # Not fixed in founders
    sum(freq > 0.03 & freq < 0.97, na.rm = TRUE) == 0,
    # Informative
    sum(freq, na.rm = TRUE) > 0.05 | sum(freq, na.rm = TRUE) < 0.95
  ) %>%
  ungroup() %>%
  # Filter to euchromatin
  filter(POS >= euchromatin_start, POS <= euchromatin_end)

cat("✓ Processed SNP data:", n_distinct(snp_data$POS), "SNPs\n")

# Apply test filters if in debug mode
if (debug_mode) {
  snp_positions <- snp_data %>%
    distinct(POS) %>%
    pull(POS)
  
  # Sample random SNPs with reproducible seed based on estimator
  # This ensures same SNPs are tested across different estimators
  seed_value <- sum(utf8ToInt(estimator)) + test_n_snps
  set.seed(seed_value)
  if (length(snp_positions) > test_n_snps) {
    test_positions <- sample(snp_positions, test_n_snps)
  } else {
    test_positions <- snp_positions
  }
  
  snp_data <- snp_data %>%
    filter(POS %in% test_positions)
  
  cat("✓ DEBUG MODE: Using", length(test_positions), "random SNPs (seed:", seed_value, ")\n")
}

# Step 2: Create base SNP dataframe with founder genotypes
cat("Creating base SNP dataframe...\n")
base_snp_df <- snp_data %>%
  # Get founder genotypes
  filter(sample_name %in% founders) %>%
  select(CHROM, POS, sample_name, freq) %>%
  pivot_wider(
    names_from = sample_name,
    values_from = freq
  ) %>%
  # Ensure founder columns are in correct order
  select(CHROM, POS, all_of(founder_order)) %>%
  # Add sample frequencies
  left_join(
    snp_data %>%
      filter(!(sample_name %in% founders)) %>%
      select(CHROM, POS, sample_name, freq) %>%
      pivot_wider(
        names_from = sample_name,
        values_from = freq
      ),
    by = c("CHROM", "POS")
  )

cat("✓ Base SNP dataframe created:", nrow(base_snp_df), "rows\n")

# Step 3: Function to interpolate haplotype frequencies for a sample
interpolate_haplotypes_for_sample <- function(haplotype_data, snp_positions, sample_name) {
  # Filter haplotype data for this sample and euchromatin
  sample_haplotypes <- haplotype_data %>%
    filter(sample == sample_name) %>%
    filter(pos >= euchromatin_start, pos <= euchromatin_end)
  
  if (nrow(sample_haplotypes) == 0) {
    return(tibble())
  }
  
  # Get haplotype positions
  haplotype_positions <- sort(unique(sample_haplotypes$pos))
  
  # For each SNP position, find interpolated haplotype frequencies
  interpolated <- tibble(
    POS = snp_positions
  ) %>%
    # Find left and right haplotype positions
    mutate(
      left_pos = map_dbl(POS, ~{
        positions <- haplotype_positions[haplotype_positions <= .x]
        if (length(positions) == 0) NA else max(positions)
      }),
      right_pos = map_dbl(POS, ~{
        positions <- haplotype_positions[haplotype_positions >= .x]
        if (length(positions) == 0) NA else min(positions)
      })
    ) %>%
    # Calculate interpolation weights
    mutate(
      alpha = (right_pos - POS) / (right_pos - left_pos),
      alpha = if_else(is.infinite(alpha) | is.nan(alpha), 0, alpha)
    ) %>%
    # Get haplotype frequencies for left and right positions
    left_join(
      sample_haplotypes %>%
        select(pos, all_of(founder_order)) %>%
        rename_with(~paste0("left_", .x), -pos) %>%
        rename(left_pos = pos),
      by = "left_pos"
    ) %>%
    left_join(
      sample_haplotypes %>%
        select(pos, all_of(founder_order)) %>%
        rename_with(~paste0("right_", .x), -pos) %>%
        rename(right_pos = pos),
      by = "right_pos"
    ) %>%
    # Interpolate each founder frequency
    mutate(across(
      starts_with("left_") & !starts_with("left_pos"),
      ~{
        founder <- str_remove(cur_column(), "left_")
        right_col <- paste0("right_", founder)
        alpha * .x + (1 - alpha) * !!sym(right_col)
      },
      .names = "{.col}"
    )) %>%
    # Select only interpolated frequencies
    select(POS, starts_with("left_") & !starts_with("left_pos")) %>%
    rename_with(~str_remove(.x, "left_"), -POS)
  
  return(interpolated)
}

# Step 4: Function to calculate imputed frequency
calculate_imputed_frequency <- function(haplotype_freqs, founder_genotypes) {
  sum(haplotype_freqs * founder_genotypes, na.rm = TRUE)
}

# Function to print compact debug output
print_compact_debug <- function(results, founder_order) {
  cat("\n=== DEBUG OUTPUT ===\n")
  
  # Header: use sprintf to control exact spacing
  header <- sprintf("%7s %32s %32s %4s %4s %3s", "POS", "Founder genotypes (%)", "Haplotype freqs (%)", "Obs", "Imp", "Dif")
  cat(header, "\n")
  
  # Founder labels: use sprintf for exact 4-character spacing
  founder_labels <- paste(sprintf("%4s", founder_order), collapse = "")
  founder_line <- sprintf("%7s %32s %32s %4s %4s %3s", "", founder_labels, founder_labels, "", "", "")
  cat(founder_line, "\n")
  
  # Dashes: use sprintf to match exact field widths
  dash_line <- sprintf("%7s %32s %32s %4s %4s %3s", 
                       strrep("-", 7), strrep("-", 32), strrep("-", 32), 
                       strrep("-", 4), strrep("-", 4), strrep("-", 3))
  cat(dash_line, "\n")
  
  for (i in 1:min(20, nrow(results))) { # Show first 20 SNPs
    row <- results[i, ]
    
    # SNP position (7 spaces)
    pos_str <- sprintf("%7d", row$POS)
    
    # Founder genotypes as percentages (4 spaces each: 3 for number + 1 for separation)
    founder_genos <- sapply(founder_order, function(f) {
      geno_val <- row[[paste0(f, "_geno")]]
      sprintf("%4.0f", round(geno_val * 100))
    })
    founder_str <- paste(founder_genos, collapse = "")
    
    # Haplotype frequencies as percentages (4 spaces each: 3 for number + 1 for separation)
    haplotype_freqs <- sapply(founder_order, function(f) {
      freq_val <- row[[paste0(f, "_freq")]]
      sprintf("%4.0f", round(freq_val * 100))
    })
    haplotype_str <- paste(haplotype_freqs, collapse = "")
    
    # Sample frequency (4 spaces, 0-100)
    sample_freq_str <- sprintf("%4.0f", round(row$observed * 100))
    
    # Imputed value (4 spaces, 0-100)
    imputed_str <- sprintf("%4.0f", round(row$imputed * 100))
    
    # Absolute difference (3 spaces, 0-100)
    diff_str <- sprintf("%3.0f", round(abs(row$observed - row$imputed) * 100))
    
              # Combine all parts with proper spacing using sprintf
          data_line <- sprintf("%7d %32s %32s %4.0f %4.0f %3.0f", 
                              row$POS, founder_str, haplotype_str, 
                              round(row$observed * 100), round(row$imputed * 100), 
                              round(abs(row$observed - row$imputed) * 100))
          cat(data_line, "\n")
  }
  
  if (nrow(results) > 20) {
    cat("... (showing first 20 of", nrow(results), "SNPs)\n")
  }
  cat("\n")
}

# Step 5: Get sample name from index
non_founder_samples <- setdiff(unique(snp_data$sample_name), founders)
if (sample_index < 1 || sample_index > length(non_founder_samples)) {
  stop("Sample index ", sample_index, " out of range. Valid range: 1-", length(non_founder_samples))
}

sample_name <- non_founder_samples[sample_index]
cat("Processing sample:", sample_name, "(index", sample_index, ")\n")

# Get interpolated haplotype frequencies for this sample
interpolated_haplotypes <- interpolate_haplotypes_for_sample(
  haplotype_results, base_snp_df$POS, sample_name
)

if (nrow(interpolated_haplotypes) == 0) {
  stop("No haplotype results for sample ", sample_name)
}

# Join with base SNP data
results <- base_snp_df %>%
  select(CHROM, POS, all_of(founder_order), !!sample_name) %>%
  rename(observed = !!sample_name) %>%
  left_join(interpolated_haplotypes, by = "POS", suffix = c("_geno", "_freq")) %>%
  # Calculate imputed frequency
  mutate(
    imputed = pmap_dbl(
      list(
        haplotype_freqs = select(cur_data(), ends_with("_freq")),
        founder_genotypes = select(cur_data(), ends_with("_geno"))
      ),
      ~calculate_imputed_frequency(as.numeric(..1), as.numeric(..2))
    )
  ) %>%
  select(CHROM, POS, all_of(paste0(founder_order, "_geno")), all_of(paste0(founder_order, "_freq")), observed, imputed) %>%
  mutate(
    sample = sample_name,
    estimator = estimator
  )

cat("  Completed:", nrow(results), "imputations\n")

# Print debug output if in debug mode
if (debug_mode) {
  print_compact_debug(results, founder_order)
}

# Step 6: Save results
if (nrow(results) > 0) {
  output_file <- if (debug_mode) {
    file.path(output_dir, paste0("snp_imputation_", estimator, "_", sample_name, "_", chr, "_DEBUG.RDS"))
  } else {
    file.path(output_dir, paste0("snp_imputation_", estimator, "_", sample_name, "_", chr, ".RDS"))
  }
  
  write_rds(results, output_file)
  
  # Summary statistics
  summary_stats <- results %>%
    summarise(
      n_imputations = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      correlation = cor(observed, imputed, use = "complete.obs"),
      rmse = sqrt(mean((observed - imputed)^2, na.rm = TRUE))
    )
  
  cat("\n✓ SNP Imputation completed successfully!\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total imputations:", nrow(results), "\n\n")
  
  cat("Summary for", sample_name, ":\n")
  print(summary_stats)
  
} else {
  cat("\n❌ No SNP imputation results obtained!\n")
}

cat("\n=== Clean SNP Imputation Complete ===\n")
