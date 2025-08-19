#!/usr/bin/env Rscript

# Test Haplotype Estimation Functions
# Professional test wrapper using the unified haplotype estimation function

source("scripts/haplotype_estimation_functions.R")

cat("=== HAPLOTYPE ESTIMATION FUNCTION TESTS ===\n\n")

# 1. Load parameters
cat("1. Loading parameters...\n")
source("helpfiles/JUICE_haplotype_parameters.R")
cat("✓ Parameter file: helpfiles/JUICE_haplotype_parameters.R\n")
cat("✓ Loaded parameters: founders (", length(founders), "), step (", step, "), h_cutoff (", h_cutoff, "), samples (", length(names_in_bam), ")\n")
cat("✓ Founders:", paste(founders, collapse=", "), "\n")

# 2. Load and transform data
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

# Transform to wide format and apply quality filter ONCE
founder_wide <- df2 %>%
  filter(name %in% founders) %>%
  select(POS, name, freq) %>%
  pivot_wider(names_from = name, values_from = freq)

# Apply quality filter: keep rows where ALL founders are fixed (< 3% or > 97%)
quality_filtered_positions <- founder_wide %>%
  filter(
    if_all(all_of(founders), ~ is.na(.x) | .x < 0.03 | .x > 0.97)
  ) %>%
  pull(POS)

# Filter to quality positions and include sample data
df3 <- df2 %>%
  filter(POS %in% quality_filtered_positions)

cat("Quality-filtered positions:", length(quality_filtered_positions), "\n")
cat("✓ Data ready:", nrow(df3), "rows\n")

# 3. Test cases
cat("\n3. Testing haplotype estimation functions...\n")

# Test parameters
test_positions <- c(5000000, 10000000, 15000000)
test_samples <- names_in_bam[1:2]  # Test first 2 samples
test_window_sizes <- c(20000, 50000)  # Fixed window sizes to test
test_h_cutoffs <- c(4, 6)  # Adaptive window h_cutoffs to test

cat("Test positions:", paste(format(test_positions, big.mark=","), collapse=", "), "\n")
cat("Test samples:", paste(test_samples, collapse=", "), "\n")
cat("Fixed window sizes:", paste(test_window_sizes/1000, "kb"), "\n")
cat("Adaptive h_cutoffs:", paste(test_h_cutoffs, collapse=", "), "\n\n")

# 4. Run tests with professional function approach
results_list <- list()

# Fixed window tests
cat("=== FIXED WINDOW TESTS ===\n")
for (window_size in test_window_sizes) {
  cat(sprintf("\n--- Testing Fixed Window: %d kb ---\n", window_size/1000))
  
  for (pos in test_positions[1:2]) {  # Test first 2 positions
    for (sample in test_samples) {
      cat(sprintf("\nTesting pos=%s, sample=%s, window=%dkb with verbose output:\n", 
                  format(pos, big.mark=","), sample, window_size/1000))
      
      # Test with detailed verbose output
      result <- estimate_haplotypes(
        pos = pos,
        sample_name = sample,
        df3 = df3,
        founders = founders,
        h_cutoff = h_cutoff,
        method = "fixed",
        window_size_bp = window_size,
        chr = "chr2R",
        verbose = 2
      )
      
      results_list[[length(results_list) + 1]] <- result
      
      cat(sprintf("Result: estimate_OK=%s, n_snps=%d\n", 
                  ifelse(is.na(result$estimate_OK), "NA", result$estimate_OK), 
                  result$n_snps))
    }
  }
}

# Adaptive window tests
cat("\n\n=== ADAPTIVE WINDOW TESTS ===\n")
for (h_cutoff_test in test_h_cutoffs) {
  cat(sprintf("\n--- Testing Adaptive Window: h_cutoff = %g ---\n", h_cutoff_test))
  
  for (pos in test_positions[1:2]) {  # Test first 2 positions
    for (sample in test_samples) {
      cat(sprintf("\nTesting pos=%s, sample=%s, h_cutoff=%g with verbose output:\n", 
                  format(pos, big.mark=","), sample, h_cutoff_test))
      
      # Test with detailed verbose output
      result <- estimate_haplotypes(
        pos = pos,
        sample_name = sample,
        df3 = df3,
        founders = founders,
        h_cutoff = h_cutoff_test,
        method = "adaptive",
        chr = "chr2R",
        verbose = 2
      )
      
      results_list[[length(results_list) + 1]] <- result
      
      cat(sprintf("Result: estimate_OK=%s, final_window=%s, n_snps=%d\n", 
                  ifelse(is.na(result$estimate_OK), "NA", result$estimate_OK),
                  ifelse(result$final_window_size >= 1000, 
                         paste0(result$final_window_size/1000, "kb"), 
                         paste0(result$final_window_size, "bp")),
                  result$n_snps))
    }
  }
}

# 5. Convert to data frame and show results
cat("\n\n=== FUNCTION TEST RESULTS SUMMARY ===\n")
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list)
  
  cat("Generated results data frame:\n")
  cat("Columns:", paste(names(results_df), collapse = ", "), "\n")
  cat("Rows:", nrow(results_df), "\n\n")
  
  # Verify all expected columns are present
  expected_cols <- c("chr", "pos", "sample", "method", "final_window_size", "n_snps", "estimate_OK", founders)
  missing_cols <- setdiff(expected_cols, names(results_df))
  extra_cols <- setdiff(names(results_df), expected_cols)
  
  if (length(missing_cols) == 0) {
    cat("✓ All expected columns present:", paste(expected_cols, collapse = ", "), "\n")
  } else {
    cat("✗ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  }
  
  if (length(extra_cols) > 0) {
    cat("ℹ Extra columns:", paste(extra_cols, collapse = ", "), "\n")
  }
  
  # Show complete results table with better formatting
  cat(sprintf("\nComplete results (%d rows):\n", nrow(results_df)))
  
  # Temporarily disable scientific notation for cleaner position display
  old_scipen <- getOption("scipen")
  options(scipen = 999)  # Avoid scientific notation
  
  print(results_df)
  
  # Restore original setting
  options(scipen = old_scipen)
  
  # Test that fixed and adaptive methods give different results
  if (any(results_df$method == "fixed") && any(results_df$method == "adaptive")) {
    cat("\n✓ Both fixed and adaptive methods tested\n")
    
    # Check that different h_cutoffs give different results for adaptive
    adaptive_results <- results_df %>% filter(method == "adaptive")
    if (length(unique(adaptive_results$h_cutoff)) > 1) {
      different_results <- adaptive_results %>%
        group_by(pos, sample) %>%
        summarise(
          n_h_cutoffs = n_distinct(h_cutoff),
          different_estimates = n_distinct(paste(estimate_OK, final_window_size)),
          .groups = "drop"
        )
      
      if (any(different_results$different_estimates > 1)) {
        cat("✓ Different h_cutoff values produce different results (constraint accumulation working)\n")
      } else {
        cat("⚠️  All h_cutoff values produce identical results (potential constraint accumulation bug)\n")
      }
    }
  }
  
  # Save test results
  test_output_file <- "test_haplotype_functions_output.RDS"
  saveRDS(results_df, test_output_file)
  cat("✓ Saved test results to:", test_output_file, "\n")
  
  # Test RDS save/load
  loaded_results <- readRDS(test_output_file)
  if (identical(results_df, loaded_results)) {
    cat("✓ RDS save/load test passed\n")
  } else {
    cat("✗ RDS save/load test failed\n")
  }
  
  # Clean up test file
  file.remove(test_output_file)
  cat("✓ Cleaned up test file\n")
  
} else {
  cat("✗ No results generated - check algorithm or data\n")
}

cat("\n=== HAPLOTYPE ESTIMATION FUNCTION TESTS COMPLETE ===\n")
