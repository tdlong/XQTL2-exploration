#!/usr/bin/env Rscript

# Quick script to peek at haplotype results
# Usage: Rscript scripts/peek_haplotype_results.R <results_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  cat("Usage: Rscript scripts/peek_haplotype_results.R <results_file>\n")
  cat("Example: Rscript scripts/peek_haplotype_results.R process/JUICE/haplotype_results/fixed_window_20kb_results_chr2R.RDS\n")
  quit(status = 1)
}

results_file <- args[1]

cat("=== PEEKING AT HAPLOTYPE RESULTS ===\n")
cat("File:", results_file, "\n\n")

# Load the results
results <- readRDS(results_file)

cat("=== BASIC INFO ===\n")
cat("Total rows:", nrow(results), "\n")
cat("Columns:", paste(names(results), collapse = ", "), "\n")
cat("Column types:", paste(sapply(results, class), collapse = ", "), "\n\n")

cat("=== SAMPLE BREAKDOWN ===\n")
sample_counts <- table(results$sample)
print(sample_counts)
cat("\n")

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

cat("=== SUCCESS RATE CALCULATION ===\n")
total_positions <- length(unique(results$pos)) * length(unique(results$sample))
cat("Total positions × samples:", total_positions, "\n")

# Count successful estimates (where B1 is not NA)
successful_estimates <- sum(!is.na(results$B1))
cat("Successful estimates (B1 not NA):", successful_estimates, "\n")

success_rate <- successful_estimates / total_positions * 100
cat("Success rate:", round(success_rate, 1), "%\n\n")

cat("=== SAMPLE-BY-SAMPLE SUCCESS ===\n")
for (sample in unique(results$sample)) {
  sample_data <- results[results$sample == sample, ]
  sample_success <- sum(!is.na(sample_data$B1))
  sample_total <- nrow(sample_data)
  sample_rate <- sample_success / sample_total * 100
  cat(sample, ":", sample_success, "/", sample_total, "(", round(sample_rate, 1), "%)\n")
}

cat("\n=== DONE ===\n")
