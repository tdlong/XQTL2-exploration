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

# Define estimators - just focus on fixed_50kb vs alternative method
estimators <- c("fixed_50kb")

# Load our haplotype estimates for all estimators
cat("Loading our haplotype estimates for all estimators...\n")

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
our_haplotypes_list <- map(estimators, load_estimator_haplotypes)
our_haplotypes_list <- our_haplotypes_list[!sapply(our_haplotypes_list, is.null)]

# Combine all estimators into one wide dataframe
cat("Combining all estimators...\n")
cat("Starting with first estimator:", names(our_haplotypes_list[[1]]), "\n")

our_haplotypes_wide <- our_haplotypes_list[[1]]
for (i in 2:length(our_haplotypes_list)) {
  cat("Joining estimator", i, ":", names(our_haplotypes_list[[i]]), "\n")
  our_haplotypes_wide <- left_join(our_haplotypes_wide, our_haplotypes_list[[i]], 
                                   by = c("chr", "pos", "estimate_OK", "sample"))
  cat("After join, columns:", names(our_haplotypes_wide), "\n")
}

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

cat("Final comparison dataframe:", nrow(comparison), "rows\n")
cat("Columns:", names(comparison), "\n\n")

# Show first 50 rows for reality check with better formatting
cat("First 50 rows of final comparison dataframe:\n")
cat("Total columns:", ncol(comparison), "\n")
cat("Column names:", paste(names(comparison), collapse = ", "), "\n\n")

# Print with options to show all columns clearly
options(width = 200)
print(head(comparison, 50), n = 50, width = Inf)

# Show summary of the two methods being compared
cat("\n=== METHOD COMPARISON SUMMARY ===\n")
cat("fixed_50kb: N =", sum(!is.na(comparison$B1_fixed_50kb)), 
    ", Mean =", round(mean(comparison$B1_fixed_50kb, na.rm = TRUE), 4),
    ", NA =", sum(is.na(comparison$B1_fixed_50kb)), "\n")
cat("B1_methodalt: N =", sum(!is.na(comparison$B1_methodalt)), 
    ", Mean =", round(mean(comparison$B1_methodalt, na.rm = TRUE), 4),
    ", NA =", sum(is.na(comparison$B1_methodalt)), "\n")

# Calculate differences between each method and alternative
cat("\n=== DIFFERENCE ANALYSIS ===\n")
diff_cols <- names(comparison)[grepl("^B1_", names(comparison)) & names(comparison) != "B1_methodalt"]

for (col in diff_cols) {
  method_name <- sub("^B1_", "", col)
  
  # Get valid pairs (both values non-NA)
  valid_pairs <- !is.na(comparison[[col]]) & !is.na(comparison$B1_methodalt)
  n_valid <- sum(valid_pairs)
  
  if (n_valid > 0) {
    # Calculate RMSE (Root Mean Square Error) only on valid pairs
    differences <- comparison[[col]][valid_pairs] - comparison$B1_methodalt[valid_pairs]
    rmse <- sqrt(mean(differences^2))
    
    # Also calculate mean absolute difference for reference
    mean_abs_diff <- mean(abs(differences))
    
    cat(sprintf("%-15s: RMSE = %.4f, Mean abs diff = %.4f (N = %d)\n", 
                method_name, rmse, mean_abs_diff, n_valid))
  } else {
    cat(sprintf("%-15s: No valid pairs for comparison (all NA)\n", method_name))
  }
}

# Show positions with large differences for any method (using RMSE-like threshold)
cat("\n=== POSITIONS WITH LARGE DIFFERENCES ===\n")
cat("Showing positions where any method differs from alternative by >0.1:\n")

large_diffs <- comparison %>%
  select(CHROM, pos, starts_with("B1_")) %>%
  mutate(across(starts_with("B1_") & !ends_with("methodalt"), 
                ~abs(. - B1_methodalt), 
                .names = "diff_{.col}")) %>%
  filter(if_any(starts_with("diff_"), ~. > 0.1)) %>%
  select(CHROM, pos, starts_with("diff_"))

if (nrow(large_diffs) > 0) {
  print(head(large_diffs, 10))
  cat("... (showing first 10 of", nrow(large_diffs), "positions with large differences)\n")
} else {
  cat("All differences are small (≤0.1) for all methods\n")
}

cat("\n=== COMPARISON COMPLETE ===\n")
