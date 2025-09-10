#!/usr/bin/env Rscript

# Analyze Combined Centromere Results
# Separates sample names and performs statistical analysis

library(tidyverse)

# Load the results
cat("Loading combined centromere results...\n")
results <- readRDS("combined_centromere_all_results.RDS")

cat("Original data structure:\n")
print(results)
cat("\nDimensions:", nrow(results), "x", ncol(results), "\n\n")

# Separate sample column into rep, TRT, sex
cat("Separating sample names...\n")
results_separated <- results %>%
  separate(sample, into = c("rep", "TRT", "sex"), sep = "_", remove = FALSE)

cat("After separation:\n")
print(head(results_separated))
cat("\n")

# Check the separation worked correctly
cat("Sample breakdown:\n")
print(table(results_separated$rep, results_separated$TRT, results_separated$sex))
cat("\n")

# Group by CHROM, pos, sex and nest
cat("Grouping by CHROM, pos, sex and nesting...\n")
nested_results <- results_separated %>%
  group_by(CHROM, pos, sex) %>%
  nest()

cat("Nested structure:\n")
print(nested_results)
cat("\n")

# Show what's in each nested group
cat("Sample of nested data (first group):\n")
if (nrow(nested_results) > 0) {
  print(nested_results$data[[1]])
}

cat("\nSummary of nested groups:\n")
nested_summary <- nested_results %>%
  mutate(n_samples = map_int(data, nrow)) %>%
  select(CHROM, pos, sex, n_samples)

print(nested_summary)

# Save the separated and nested results
saveRDS(results_separated, "combined_centromere_results_separated.RDS")
saveRDS(nested_results, "combined_centromere_results_nested.RDS")

cat("\n✓ Results saved:\n")
cat("- combined_centromere_results_separated.RDS\n")
cat("- combined_centromere_results_nested.RDS\n")
