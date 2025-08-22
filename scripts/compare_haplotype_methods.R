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

# Define all estimators
estimators <- c(
  "fixed_20kb", "fixed_50kb", "fixed_100kb", "fixed_200kb", "fixed_500kb",
  "adaptive_h4", "adaptive_h6", "adaptive_h8", "adaptive_h10"
)

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
  
  if (length(target_index) == 0) {
    stop("Sample", target_sample, "not found in alternative data")
  }
  
  # Extract B1 frequencies for our target sample
  B1_freqs <- map_dbl(1:nrow(alt_data), function(i) {
    # Get haplotype frequencies for this position and sample
    hap_freqs <- alt_data$Haps[[i]][[1]][[target_index]]
    # Return B1 frequency
    hap_freqs["B1"]
  })
  
  # Create clean dataframe
  result <- alt_data %>%
    select(CHROM, pos) %>%
    mutate(B1_methodalt = B1_freqs)
  
  return(result)
}

alt_B1_freqs <- extract_B1_freqs(alt_haplotypes, sample_name)
cat("Extracted B1 frequencies for", nrow(alt_B1_freqs), "positions\n\n")

# Join our estimators with alternative method
cat("=== COMBINED COMPARISON ===\n")
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

# Print with options to show all columns
options(width = 200)
print(head(comparison, 50), n = 50, width = Inf)

# Also show summary of each column
cat("\n=== COLUMN SUMMARIES ===\n")
for (col in names(comparison)) {
  if (grepl("^B1_", col)) {
    cat(sprintf("%-20s: N = %d, Mean = %.4f, NA = %d\n", 
                col, 
                sum(!is.na(comparison[[col]])), 
                mean(comparison[[col]], na.rm = TRUE),
                sum(is.na(comparison[[col]]))))
  }
}

# Calculate differences between each method and alternative
cat("\n=== DIFFERENCE ANALYSIS ===\n")
diff_cols <- names(comparison)[grepl("^B1_", names(comparison)) & names(comparison) != "B1_methodalt"]

for (col in diff_cols) {
  method_name <- sub("^B1_", "", col)
  
  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((comparison[[col]] - comparison$B1_methodalt)^2, na.rm = TRUE))
  
  # Also calculate mean absolute difference for reference
  mean_abs_diff <- mean(abs(comparison[[col]] - comparison$B1_methodalt), na.rm = TRUE)
  
  cat(sprintf("%-15s: RMSE = %.4f, Mean abs diff = %.4f\n", method_name, rmse, mean_abs_diff))
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
