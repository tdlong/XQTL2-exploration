#!/usr/bin/env Rscript

# Debug script to check what haplotype files exist
cat("=== CHECKING HAPLOTYPE FILES ===\n")

# Check our haplotype results directory
cat("Checking process/JUICE/haplotype_results/:\n")
if (dir.exists("process/JUICE/haplotype_results/")) {
  files <- list.files("process/JUICE/haplotype_results/", pattern = "*_results_chr2R.RDS")
  cat("Haplotype files found:", length(files), "\n")
  print(files)
} else {
  cat("Directory does not exist\n")
}

cat("\n=== CHECKING ALTERNATIVE FILE ===\n")
alt_file <- "/dfs7/adl/tdlong/fly_pool/XQTL2/process/JUICE/R.haps.chr2R.out.rds"
if (file.exists(alt_file)) {
  cat("Alternative file exists:", alt_file, "\n")
  cat("File size:", file.size(alt_file), "bytes\n")
} else {
  cat("Alternative file does not exist:", alt_file, "\n")
}

cat("\n=== CHECKING CURRENT WORKING DIRECTORY ===\n")
cat("Current working directory:", getwd(), "\n")
cat("Files in current directory:\n")
print(list.files(".", pattern = "*.RDS"))

