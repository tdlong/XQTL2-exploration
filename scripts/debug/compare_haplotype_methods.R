#!/usr/bin/env Rscript

# Compare haplotype frequencies between all our estimators and alternative method
# Extract B1 frequencies for AJ_1_1 from all methods

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Load parameter file to get sample names
source("helpfiles/JUICE_haplotype_parameters.R")
sample_name <- "AJ_1_1"
sample_index <- which(names_in_bam == sample_name)

cat("=== Comparing All Haplotype Methods ===\n")
cat("Sample:", sample_name, "(index:", sample_index, ")\n\n")

# Load all estimators to get the full comparison dataframe
cat("Loading our haplotype estimates for all estimators...\n")

# Define all estimators
all_estimators <- c(
  "fixed_20kb", "fixed_50kb", "fixed_100kb", "fixed_200kb", "fixed_500kb",
  "adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10"
)

load_estimator_haplotypes <- function(estimator) {
  if (grepl("^fixed_", estimator)) {
    window_size <- as.numeric(sub("fixed_(\\d+)kb", "\\1", estimator))
    file_path <- paste0("process/JUICE/haplotype_results/fixed_window_", window_size, "kb_results_chr2R.RDS")
  } else if (grepl("^adaptive_", estimator)) {
    h_cutoff <- as.numeric(sub("adaptive_h(\\d+)", "\\1", estimator))
    file_path <- paste0("process/JUICE/haplotype_results/adaptive_window_h", h_cutoff, "_results_chr2R.RDS")
  }
  
  if (file.exists(file_path)) {
    # Copy exact working code from create_summary_file_chunked.R
    data <- readRDS(file_path) %>%
      select(chr, pos, estimate_OK, B1, sample) %>%
      filter(sample == sample_name)
    
    # Rename B1 to include estimator name
    names(data)[names(data) == "B1"] <- paste0("B1_", estimator)
    return(data)
  } else {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
}

# Load all estimators
our_haplotypes_list <- map(all_estimators, load_estimator_haplotypes)
our_haplotypes_list <- our_haplotypes_list[!sapply(our_haplotypes_list, is.null)]

# Combine all estimators into one wide dataframe
cat("Combining all estimators...\n")
our_haplotypes_wide <- reduce(our_haplotypes_list, function(x, y) {
  left_join(x, y, by = c("chr", "pos", "estimate_OK", "sample"))
}, .init = our_haplotypes_list[[1]])

cat("Our estimates combined:", nrow(our_haplotypes_wide), "positions\n")
cat("Columns:", names(our_haplotypes_wide), "\n\n")

# Load alternative haplotype estimates
cat("Loading alternative haplotype estimates...\n")
alt_haplotypes <- readRDS("/dfs7/adl/tdlong/fly_pool/XQTL2/process/JUICE/R.haps.chr2R.out.rds")

cat("Alternative estimates:", nrow(alt_haplotypes), "positions\n")
cat("Sample names in alternative data:", unlist(alt_haplotypes$sample[[1]]), "\n\n")

# Extract B1 frequencies for AJ_1_1 from the nested structure
cat("Extracting B1 frequencies for", sample_name, "from alternative method...\n")

extract_B1_freqs <- function(alt_data, target_sample) {
  # Find the index of our target sample
  sample_names <- unlist(alt_data$sample[[1]])
  target_index <- which(sample_names == target_sample)
  
  cat("Target sample:", target_sample, "found at index:", target_index, "\n")
  cat("Sample names:", paste(sample_names, collapse = ", "), "\n")
  
  if (length(target_index) == 0) {
    stop("Sample", target_sample, "not found in alternative data")
  }
  
  # Debug the first few extractions
  cat("Debugging first 5 extractions:\n")
  for (i in 1:5) {
    cat("  Position", i, ":", alt_data$pos[i], "\n")
    cat("    Haps structure:", str(alt_data$Haps[[i]]), "\n")
    if (length(alt_data$Haps[[i]]) > 0) {
      # Try to access the haplotype data for our target sample
      hap_freqs <- alt_data$Haps[[i]][[target_index]]
      b1_freq <- hap_freqs[1]  # First founder frequency (B1)
      cat("    -> B1 =", b1_freq, "\n")
      cat("    -> B1 class:", class(b1_freq), "\n")
      cat("    -> Full haplotype vector:", paste(names(hap_freqs), "=", hap_freqs, collapse = ", "), "\n")
    } else {
      cat("    -> No haplotype data\n")
    }
  }
  
  # Extract B1 frequencies for our target sample
  B1_freqs <- map_dbl(1:nrow(alt_data), function(i) {
    # Get haplotype frequencies for this position and sample
    if (length(alt_data$Haps[[i]]) > 0) {
      hap_freqs <- alt_data$Haps[[i]][[target_index]]
      # Extract B1 (first founder frequency) and convert to numeric
      b1_value <- hap_freqs[1]
      return(as.numeric(b1_value))
    } else {
      return(NA_real_)
    }
  })
  
  cat("Extracted B1 frequencies summary:\n")
  cat("  Total positions:", length(B1_freqs), "\n")
  cat("  Non-NA values:", sum(!is.na(B1_freqs)), "\n")
  cat("  NA values:", sum(is.na(B1_freqs)), "\n")
  cat("  Range:", range(B1_freqs, na.rm = TRUE), "\n")
  
  # Create clean dataframe
  result <- alt_data %>%
    select(CHROM, pos) %>%
    mutate(B1_methodalt = B1_freqs)
  
  cat("Result dataframe created:", nrow(result), "rows\n")
  cat("Sample result rows:\n")
  print(head(result, 5))
  
  return(result)
}

alt_B1_freqs <- extract_B1_freqs(alt_haplotypes, sample_name)
cat("Extracted B1 frequencies for", nrow(alt_B1_freqs), "positions\n\n")

# Debug position ranges before join
cat("=== POSITION RANGE DEBUG ===\n")
cat("Our haplotype positions range:", range(our_haplotypes_wide$pos), "\n")
cat("Alternative method positions range:", range(alt_B1_freqs$pos), "\n")
cat("Position overlap check:\n")
cat("  Our min pos:", min(our_haplotypes_wide$pos), "vs Alt min pos:", min(alt_B1_freqs$pos), "\n")
cat("  Our max pos:", max(our_haplotypes_wide$pos), "vs Alt max pos:", max(alt_B1_freqs$pos), "\n")

# Check for exact position matches
common_positions <- intersect(our_haplotypes_wide$pos, alt_B1_freqs$pos)
cat("Common positions:", length(common_positions), "out of", nrow(our_haplotypes_wide), "\n")

# Show sample positions from each dataset
cat("Sample positions from our data:", head(our_haplotypes_wide$pos, 10), "\n")
cat("Sample positions from alt data:", head(alt_B1_freqs$pos, 10), "\n")

# Join our estimators with alternative method
cat("\n=== COMBINED COMPARISON ===\n")
comparison <- our_haplotypes_wide %>%
  left_join(alt_B1_freqs, by = "pos") %>%
  select(CHROM, pos, starts_with("B1_")) %>%
  arrange(pos)

# Now subset to just the three methods we want to compare
comparison_subset <- comparison %>%
  select(CHROM, pos, B1_fixed_50kb, B1_fixed_100kb, B1_methodalt) %>%
  filter(!is.na(B1_fixed_50kb) & !is.na(B1_fixed_100kb) & !is.na(B1_methodalt))

cat("Subset for comparison (fixed_50kb, fixed_100kb vs alternative):", nrow(comparison_subset), "positions\n")

cat("Final comparison dataframe:", nrow(comparison), "rows\n")
cat("Columns:", names(comparison), "\n\n")

# Show the subset for the three methods we want to compare
cat("First 50 rows of comparison subset (fixed_50kb, fixed_100kb vs alternative):\n")
cat("Columns:", paste(names(comparison_subset), collapse = ", "), "\n\n")

# Print the subset clearly
options(width = 200)
print(head(comparison_subset, 50), n = 50, width = Inf)

# Show summary of the three methods being compared
cat("\n=== METHOD COMPARISON SUMMARY ===\n")
cat("fixed_50kb: N =", sum(!is.na(comparison_subset$B1_fixed_50kb)), 
    ", Mean =", round(mean(comparison_subset$B1_fixed_50kb, na.rm = TRUE), 4),
    ", NA =", sum(is.na(comparison_subset$B1_fixed_50kb)), "\n")
cat("fixed_100kb: N =", sum(!is.na(comparison_subset$B1_fixed_100kb)), 
    ", Mean =", round(mean(comparison_subset$B1_fixed_100kb, na.rm = TRUE), 4),
    ", NA =", sum(is.na(comparison_subset$B1_fixed_100kb)), "\n")
cat("B1_methodalt: N =", sum(!is.na(comparison_subset$B1_methodalt)), 
    ", Mean =", round(mean(comparison_subset$B1_methodalt, na.rm = TRUE), 4),
    ", NA =", sum(is.na(comparison_subset$B1_methodalt)), "\n")

# Calculate differences between both estimators and alternative method
cat("\n=== DIFFERENCE ANALYSIS ===\n")
cat("Comparing both estimators vs alternative method:\n")

# Calculate RMSE and mean absolute difference for fixed_50kb
diff_50kb <- comparison_subset$B1_fixed_50kb - comparison_subset$B1_methodalt
rmse_50kb <- sqrt(mean(diff_50kb^2))
mean_abs_diff_50kb <- mean(abs(diff_50kb))

# Calculate RMSE and mean absolute difference for fixed_100kb
diff_100kb <- comparison_subset$B1_fixed_100kb - comparison_subset$B1_methodalt
rmse_100kb <- sqrt(mean(diff_100kb^2))
mean_abs_diff_100kb <- mean(abs(diff_100kb))

cat(sprintf("fixed_50kb vs alternative: RMSE = %.4f, Mean abs diff = %.4f (N = %d)\n", 
            rmse_50kb, mean_abs_diff_50kb, nrow(comparison_subset)))
cat(sprintf("fixed_100kb vs alternative: RMSE = %.4f, Mean abs diff = %.4f (N = %d)\n", 
            rmse_100kb, mean_abs_diff_100kb, nrow(comparison_subset)))

# Show positions with large differences between both estimators and alternative
cat("\n=== POSITIONS WITH LARGE DIFFERENCES ===\n")
cat("Showing positions where either estimator differs from alternative by >0.1:\n")

# Check for large differences for both estimators
large_diff_positions <- comparison_subset %>%
  mutate(
    diff_50kb = abs(B1_fixed_50kb - B1_methodalt),
    diff_100kb = abs(B1_fixed_100kb - B1_methodalt)
  ) %>%
  filter(diff_50kb > 0.1 | diff_100kb > 0.1) %>%
  arrange(desc(pmax(diff_50kb, diff_100kb)))

if (nrow(large_diff_positions) > 0) {
  cat("Found", nrow(large_diff_positions), "positions with large differences:\n")
  print(large_diff_positions, n = min(20, nrow(large_diff_positions)))
} else {
  cat("All differences are small (≤0.1) between both estimators and alternative\n")
}

cat("\n=== COMPARISON COMPLETE ===\n")
