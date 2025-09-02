#!/usr/bin/env Rscript

# Clean SNP Imputation - Tidyverse Approach
# Based on clean dataframe structure with founder genotypes and haplotype frequencies
# 
# Usage: Rscript snp_imputation_clean.R <chr> <param_file> <output_dir> <estimator> <sample_index> [test_n_snps]

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6) {
  stop("Usage: Rscript snp_imputation_clean.R <chr> <param_file> <output_dir> <estimator> <sample_index> [test_n_snps]")
}

chr <- args[1]
param_file <- args[2]
output_dir <- args[3]
estimator <- args[4]
sample_index <- as.numeric(args[5])
test_n_snps <- if (length(args) >= 6) as.numeric(args[6]) else NULL

# Determine if we're in debug mode
debug_mode <- !is.null(test_n_snps)

cat("=== Clean SNP Imputation ===\n")
cat("Chromosome:", chr, "\n")
cat("Estimator:", estimator, "\n")
cat("Sample index:", sample_index, "\n")
if (debug_mode) {
  cat("DEBUG MODE: Testing with", test_n_snps, "random SNPs\n")
}
cat("\n")

# Load parameter file
cat("Loading parameter file...\n")
source(param_file)
cat("✓ Parameter file loaded\n")
cat("Founders:", paste(founders, collapse = ", "), "\n")

# Get sample name from index
if (sample_index < 1 || sample_index > length(names_in_bam)) {
  stop("Invalid sample index. Must be between 1 and", length(names_in_bam))
}
sample_name <- names_in_bam[sample_index]
cat("Sample name:", sample_name, "\n\n")

# Define euchromatin boundaries for each chromosome
euchromatin_boundaries <- list(
  chr2L = c(82455, 22011009),
  chr2R = c(5398184, 24684540),
  chr3L = c(158639, 22962476),
  chr3R = c(4552934, 31845060),
  chrX = c(277911, 22628490)
)

# Get boundaries for this chromosome
if (!chr %in% names(euchromatin_boundaries)) {
  stop("Invalid chromosome: ", chr, ". Valid chromosomes: ", paste(names(euchromatin_boundaries), collapse = ", "))
}

euchromatin_start <- euchromatin_boundaries[[chr]][1]
euchromatin_end <- euchromatin_boundaries[[chr]][2]
cat("Euchromatin region for", chr, ":", euchromatin_start, "-", euchromatin_end, "bp\n\n")

# Load haplotype results for this estimator
cat("Loading haplotype estimation results...\n")
if (grepl("^fixed_", estimator)) {
  # Fixed window estimator
  window_size_kb <- as.numeric(gsub("fixed_", "", gsub("kb", "", estimator)))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("fixed_window_", window_size_kb, "kb_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file)
  cat("✓ Fixed window results loaded for", window_size_kb, "kb window\n")
} else if (grepl("^adaptive_h", estimator)) {
  # Adaptive window estimator
  h_cutoff <- as.numeric(gsub("adaptive_h", "", estimator))
  haplotype_file <- file.path(output_dir, "haplotype_results", paste0("adaptive_window_h", h_cutoff, "_results_", chr, ".RDS"))
  haplotype_results <- read_rds(haplotype_file)
  cat("✓ Adaptive window results loaded for h_cutoff =", h_cutoff, "\n")
} else {
  stop("Invalid estimator format. Must be 'fixed_<size>kb' or 'adaptive_h<number>'")
}

cat("✓ Haplotype results loaded:", nrow(haplotype_results), "rows\n")

# Load observed SNP data from REFALT files (keep working preprocessing)
cat("Loading observed SNP data from REFALT files...\n")
refalt_file <- file.path(output_dir, paste0("RefAlt.", chr, ".txt"))
if (!file.exists(refalt_file)) {
  stop("REFALT data file not found: ", refalt_file)
}

# Load and process REFALT data (keep working preprocessing)
df <- read.table(refalt_file, header = TRUE)
cat("✓ Raw REFALT data loaded:", nrow(df), "rows\n")

# Transform REF/ALT counts to frequencies (keep working preprocessing)
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

cat("✓ Processed REFALT data:", nrow(df2), "rows\n")

# Filter for high-quality SNPs (keep working preprocessing)
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

# Get valid SNPs for evaluation (euchromatin only)
valid_snps <- good_snps %>%
  filter(POS >= euchromatin_start & POS <= euchromatin_end) %>%
  left_join(df2, multiple = "all")

cat("✓ Valid euchromatic SNPs for evaluation:", nrow(valid_snps %>% distinct(CHROM, POS)), "\n\n")

# Apply test filters if in debug mode
if (debug_mode) {
  snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
  
  # Sample random SNPs with reproducible seed based on estimator
  seed_value <- sum(utf8ToInt(estimator)) + test_n_snps
  set.seed(seed_value)
  if (length(snp_positions) > test_n_snps) {
    test_positions <- sample(snp_positions, test_n_snps)
  } else {
    test_positions <- snp_positions
  }
  
  valid_snps <- valid_snps %>% filter(POS %in% test_positions)
  
  cat("✓ DEBUG MODE: Using", length(test_positions), "random SNPs (seed:", seed_value, ")\n")
}

# Get SNP positions for this sample
snp_positions <- valid_snps %>% distinct(CHROM, POS) %>% pull(POS)
cat("SNP positions:", length(snp_positions), "\n\n")

# Get haplotype results for this sample
sample_haplotypes <- haplotype_results %>%
  filter(sample == sample_name)

if (nrow(sample_haplotypes) == 0) {
  stop("No haplotype results for sample", sample_name)
}

# Filter haplotype results to euchromatin
sample_haplotypes <- sample_haplotypes %>%
  filter(pos >= euchromatin_start & pos <= euchromatin_end)

cat("Haplotype positions for sample:", nrow(sample_haplotypes), "\n")

# Define founder order consistently
founder_order <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")

# Create clean dataframe structure with founder genotypes and haplotype frequencies
cat("Creating clean dataframe structure...\n")

# Step 1: Get founder genotypes for each SNP position
founder_genotypes <- valid_snps %>%
  filter(name %in% founders) %>%
  select(CHROM, POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq) %>%
  # Ensure founder columns are in correct order
  select(CHROM, POS, all_of(founder_order)) %>%
  # Rename columns to indicate these are genotypes
  rename_with(~paste0(.x, "_geno"), -c(CHROM, POS))

# Step 2: Get observed frequencies for this sample
sample_observed <- valid_snps %>%
  filter(name == sample_name) %>%
  select(CHROM, POS, freq) %>%
  rename(observed = freq)

# Step 3: Interpolate haplotype frequencies to SNP positions
cat("Interpolating haplotype frequencies to SNP positions...\n")

# Get haplotype positions (sorted)
haplotype_positions <- sort(unique(sample_haplotypes$pos))

# Create interpolation function
interpolate_haplotypes <- function(snp_pos, haplotype_data, founder_order) {
  # Find left and right haplotype positions
  left_pos <- haplotype_positions[haplotype_positions <= snp_pos]
  right_pos <- haplotype_positions[haplotype_positions >= snp_pos]
  
  if (length(left_pos) == 0 || length(right_pos) == 0) {
    return(rep(NA, length(founder_order)))
  }
  
  left_pos <- max(left_pos)
  right_pos <- min(right_pos)
  
  # If SNP is exactly at a haplotype position, use that
  if (left_pos == right_pos) {
    haplotype_row <- haplotype_data %>% filter(pos == left_pos)
    if (nrow(haplotype_row) == 0) return(rep(NA, length(founder_order)))
    return(as.numeric(haplotype_row[1, founder_order]))
  }
  
  # Linear interpolation
  left_freqs <- haplotype_data %>% filter(pos == left_pos) %>% select(all_of(founder_order))
  right_freqs <- haplotype_data %>% filter(pos == right_pos) %>% select(all_of(founder_order))
  
  if (nrow(left_freqs) == 0 || nrow(right_freqs) == 0) {
    return(rep(NA, length(founder_order)))
  }
  
  left_freqs <- as.numeric(left_freqs[1, ])
  right_freqs <- as.numeric(right_freqs[1, ])
  
  # Interpolation weight
  alpha <- (right_pos - snp_pos) / (right_pos - left_pos)
  
  # Interpolate each founder frequency
  interpolated_freqs <- alpha * left_freqs + (1 - alpha) * right_freqs
  
  return(interpolated_freqs)
}

# Apply interpolation to each SNP position
interpolated_haplotypes <- map_dfr(snp_positions, function(pos) {
  freqs <- interpolate_haplotypes(pos, sample_haplotypes, founder_order)
  tibble(
    POS = pos,
    !!!setNames(as.list(freqs), paste0(founder_order, "_freq"))
  )
}, .id = NULL)

# Step 4: Combine everything into clean dataframe
cat("Combining into clean dataframe structure...\n")

results <- founder_genotypes %>%
  left_join(sample_observed, by = c("CHROM", "POS")) %>%
  left_join(interpolated_haplotypes, by = "POS") %>%
  # Calculate imputed frequency: sum(haplotype_freqs * founder_genotypes)
  rowwise() %>%
  mutate(
    imputed = sum(
      B1_freq * B1_geno + B2_freq * B2_geno + B3_freq * B3_geno + B4_freq * B4_geno +
      B5_freq * B5_geno + B6_freq * B6_geno + B7_freq * B7_geno + AB8_freq * AB8_geno,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  # Add sample and estimator info
  mutate(
    sample = sample_name,
    estimator = estimator
  ) %>%
  # Select columns in logical order
  select(CHROM, POS, sample, estimator, observed, imputed, 
         ends_with("_geno"), ends_with("_freq"))

cat("✓ Clean dataframe created:", nrow(results), "rows\n")

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

# Print debug output if in debug mode
if (debug_mode && nrow(results) > 0) {
  print_compact_debug(results, founder_order)
}

# Save results
if (nrow(results) > 0) {
  output_file <- if (debug_mode) {
    file.path(output_dir, "haplotype_results", paste0("snp_imputation_", estimator, "_", sample_name, "_", chr, "_DEBUG.RDS"))
  } else {
    file.path(output_dir, "haplotype_results", paste0("snp_imputation_", estimator, "_", sample_name, "_", chr, ".RDS"))
  }
  
  write_rds(results, output_file)
  
  # Summary statistics
  summary_stats <- results %>%
    summarize(
      n_imputations = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_imputed = mean(imputed, na.rm = TRUE),
      correlation = cor(observed, imputed, use = "complete.obs"),
      rmse = sqrt(mean((observed - imputed)^2, na.rm = TRUE))
    )
  
  cat("\n✓ SNP Imputation completed successfully!\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total imputations:", nrow(results), "\n")
  cat("Correlation:", round(summary_stats$correlation, 4), "\n")
  cat("RMSE:", round(summary_stats$rmse, 4), "\n")
  
} else {
  cat("\n❌ No SNP imputation results obtained!\n")
}

cat("\n=== SNP Imputation Complete ===\n")
