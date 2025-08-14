#!/usr/bin/env Rscript

# Simple test script to run haplotype estimation for chr2R fixed 20kb
# and then peek at the results

cat("=== TESTING SINGLE HAPLOTYPE ESTIMATION ===\n")
cat("Chromosome: chr2R\n")
cat("Method: fixed\n")
cat("Parameter: 20kb\n\n")

# Set parameters
chr <- "chr2R"
method <- "fixed"
param <- 20
param_file <- "helpfiles/JUICE/JUICE_haplotype_parameters.R"
output_dir <- "process/JUICE"

cat("Running haplotype estimation...\n")

# Run the haplotype estimation
if (method == "fixed") {
  system(paste("Rscript scripts/REFALT2haps.FixedWindow.Single.R", chr, param_file, output_dir, param))
} else {
  system(paste("Rscript scripts/REFALT2haps.AdaptWindow.Single.R", chr, param_file, output_dir, param))
}

cat("\n=== HAPLOTYPE ESTIMATION COMPLETE ===\n")

# Now peek at the results
results_file <- paste0(output_dir, "/haplotype_results/fixed_window_", param, "kb_results_", chr, ".RDS")

if (file.exists(results_file)) {
  cat("Results file found:", results_file, "\n\n")
  
  # Load and examine results
  results <- readRDS(results_file)
  
  cat("=== BASIC INFO ===\n")
  cat("Total rows:", nrow(results), "\n")
  cat("Columns:", paste(names(results), collapse = ", "), "\n")
  cat("Column types:", paste(sapply(results, class), collapse = ", "), "\n\n")
  
  cat("=== FIRST FEW ROWS ===\n")
  print(head(results, 3))
  cat("\n")
  
  cat("=== FOUNDER FREQUENCY STATS ===\n")
  founder_cols <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")
  founder_cols <- founder_cols[founder_cols %in% names(results)]
  
  for (col in founder_cols) {
    cat(col, ":\n")
    cat("  Non-NA values:", sum(!is.na(results[[col]])), "\n")
    cat("  NA values:", sum(is.na(results[[col]])), "\n")
    if (sum(!is.na(results[[col]])) > 0) {
      cat("  Range:", range(results[[col]], na.rm = TRUE), "\n")
      cat("  Mean:", mean(results[[col]], na.rm = TRUE), "\n")
    }
    cat("\n")
  }
  
  cat("=== SUCCESS RATE ===\n")
  total_positions <- length(unique(results$pos)) * length(unique(results$sample))
  successful_estimates <- sum(!is.na(results$B1))
  success_rate <- successful_estimates / total_positions * 100
  
  cat("Total positions × samples:", total_positions, "\n")
  cat("Successful estimates (B1 not NA):", successful_estimates, "\n")
  cat("Success rate:", round(success_rate, 1), "%\n\n")
  
} else {
  cat("❌ Results file not found:", results_file, "\n")
}

cat("=== TEST COMPLETE ===\n")
