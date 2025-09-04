#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

# =============================================================================
# Check Adaptive H4 List Format Results
# =============================================================================
# 
# This script checks the adaptive_h4 list format results to see what's actually
# being saved and diagnose why all positions might have groups = rep(1,8),
# haps = NA, and error matrices = NA.
#
# USAGE:
# Rscript scripts/debug/check_adaptive_h4_list_results.R <chr> <output_dir>
#
# EXAMPLE:
# Rscript scripts/debug/check_adaptive_h4_list_results.R chr2R process/ZINC2
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/debug/check_adaptive_h4_list_results.R <chr> <output_dir>")
}

chr <- args[1]
output_dir <- args[2]

cat("=== CHECKING ADAPTIVE H4 LIST FORMAT RESULTS ===\n")
cat("Chromosome:", chr, "\n")
cat("Output directory:", output_dir, "\n\n")

# Check if the file exists (try both old and new directory structures)
adaptive_file_new <- file.path(output_dir, "list_results", paste0("adaptive_h4_list_format_", chr, ".RDS"))
adaptive_file_old <- file.path(output_dir, "hap_list_results", paste0("adaptive_h4_list_format_", chr, ".RDS"))

if (file.exists(adaptive_file_new)) {
  adaptive_file <- adaptive_file_new
  cat("Using new directory structure: list_results/\n")
} else if (file.exists(adaptive_file_old)) {
  adaptive_file <- adaptive_file_old
  cat("Using old directory structure: hap_list_results/\n")
} else {
  cat("❌ Adaptive h4 results file not found in either location:\n")
  cat("  New:", adaptive_file_new, "\n")
  cat("  Old:", adaptive_file_old, "\n")
  quit(status = 1)
}


cat("✓ Loading adaptive_h4 results from:", adaptive_file, "\n")
adaptive_data <- readRDS(adaptive_file)

cat("✓ Data loaded successfully\n")
cat("Dimensions:", nrow(adaptive_data), "rows,", ncol(adaptive_data), "columns\n")
cat("Columns:", paste(names(adaptive_data), collapse = ", "), "\n\n")

# Check the structure
cat("=== DATA STRUCTURE ===\n")
str(adaptive_data, max.level = 2)

# Check first few rows
cat("\n=== FIRST 5 ROWS ===\n")
print(head(adaptive_data, 5))

# Check groups
cat("\n=== GROUPS ANALYSIS ===\n")
if ("Groups" %in% names(adaptive_data)) {
  # Check first few groups
  cat("First 3 Groups entries:\n")
  for (i in 1:min(3, nrow(adaptive_data))) {
    cat("Position", adaptive_data$pos[i], ":\n")
    groups_list <- adaptive_data$Groups[[i]]
    if (!is.null(groups_list) && length(groups_list) > 0) {
      for (j in 1:min(2, length(groups_list))) {
        groups <- groups_list[[j]]
        cat("  Sample", j, "groups:", paste(groups, collapse = ", "), "\n")
        cat("  Unique groups:", length(unique(groups)), "\n")
        cat("  Is 1:8?", all(sort(unique(groups)) == 1:8), "\n")
      }
    } else {
      cat("  Groups is NULL or empty\n")
    }
  }
  
  # Count how many positions have groups = 1:8 vs other
  groups_analysis <- adaptive_data %>%
    mutate(
      has_groups = map_lgl(Groups, ~ !is.null(.x) && length(.x) > 0),
      first_sample_groups = map(Groups, ~ if (!is.null(.x) && length(.x) > 0) .x[[1]] else NA),
      is_1_to_8 = map_lgl(first_sample_groups, ~ {
        if (is.na(.x) || is.null(.x)) return(FALSE)
        all(sort(unique(.x)) == 1:8)
      }),
      n_unique_groups = map_dbl(first_sample_groups, ~ {
        if (is.na(.x) || is.null(.x)) return(NA)
        length(unique(.x))
      })
    )
  
  cat("\nGroups summary:\n")
  cat("Positions with groups data:", sum(groups_analysis$has_groups, na.rm = TRUE), "/", nrow(groups_analysis), "\n")
  cat("Positions with groups = 1:8:", sum(groups_analysis$is_1_to_8, na.rm = TRUE), "/", nrow(groups_analysis), "\n")
  
  if (sum(groups_analysis$has_groups, na.rm = TRUE) > 0) {
    cat("Unique groups distribution:\n")
    print(table(groups_analysis$n_unique_groups, useNA = "ifany"))
  }
}

# Check haplotypes
cat("\n=== HAPLOTYPES ANALYSIS ===\n")
if ("Haps" %in% names(adaptive_data)) {
  # Check first few haplotypes
  cat("First 3 Haps entries:\n")
  for (i in 1:min(3, nrow(adaptive_data))) {
    cat("Position", adaptive_data$pos[i], ":\n")
    haps_list <- adaptive_data$Haps[[i]]
    if (!is.null(haps_list) && length(haps_list) > 0) {
      for (j in 1:min(2, length(haps_list))) {
        haps <- haps_list[[j]]
        if (!is.null(haps)) {
          cat("  Sample", j, "haps:", paste(round(haps, 3), collapse = ", "), "\n")
          cat("  Sum:", round(sum(haps, na.rm = TRUE), 3), "\n")
          cat("  All NA?", all(is.na(haps)), "\n")
        } else {
          cat("  Sample", j, "haps is NULL\n")
        }
      }
    } else {
      cat("  Haps is NULL or empty\n")
    }
  }
  
  # Count how many positions have NA haplotypes
  haps_analysis <- adaptive_data %>%
    mutate(
      has_haps = map_lgl(Haps, ~ !is.null(.x) && length(.x) > 0),
      first_sample_haps = map(Haps, ~ if (!is.null(.x) && length(.x) > 0) .x[[1]] else NA),
      all_na_haps = map_lgl(first_sample_haps, ~ {
        if (is.na(.x) || is.null(.x)) return(TRUE)
        all(is.na(.x))
      })
    )
  
  cat("\nHaplotypes summary:\n")
  cat("Positions with haps data:", sum(haps_analysis$has_haps, na.rm = TRUE), "/", nrow(haps_analysis), "\n")
  cat("Positions with all NA haps:", sum(haps_analysis$all_na_haps, na.rm = TRUE), "/", nrow(haps_analysis), "\n")
}

# Check error matrices
cat("\n=== ERROR MATRICES ANALYSIS ===\n")
if ("Err" %in% names(adaptive_data)) {
  # Check first few error matrices
  cat("First 3 Err entries:\n")
  for (i in 1:min(3, nrow(adaptive_data))) {
    cat("Position", adaptive_data$pos[i], ":\n")
    err_list <- adaptive_data$Err[[i]]
    if (!is.null(err_list) && length(err_list) > 0) {
      for (j in 1:min(2, length(err_list))) {
        err <- err_list[[j]]
        if (!is.null(err)) {
          cat("  Sample", j, "error matrix dimensions:", nrow(err), "x", ncol(err), "\n")
          cat("  All NA?", all(is.na(err)), "\n")
          if (!all(is.na(err))) {
            cat("  Non-NA values:", sum(!is.na(err)), "/", length(err), "\n")
          }
        } else {
          cat("  Sample", j, "error matrix is NULL\n")
        }
      }
    } else {
      cat("  Err is NULL or empty\n")
    }
  }
  
  # Count how many positions have NA error matrices
  err_analysis <- adaptive_data %>%
    mutate(
      has_err = map_lgl(Err, ~ !is.null(.x) && length(.x) > 0),
      first_sample_err = map(Err, ~ if (!is.null(.x) && length(.x) > 0) .x[[1]] else NA),
      all_na_err = map_lgl(first_sample_err, ~ {
        if (is.na(.x) || is.null(.x)) return(TRUE)
        all(is.na(.x))
      })
    )
  
  cat("\nError matrices summary:\n")
  cat("Positions with err data:", sum(err_analysis$has_err, na.rm = TRUE), "/", nrow(err_analysis), "\n")
  cat("Positions with all NA err:", sum(err_analysis$all_na_err, na.rm = TRUE), "/", nrow(err_analysis), "\n")
}

# Check samples
cat("\n=== SAMPLES ANALYSIS ===\n")
if ("sample" %in% names(adaptive_data)) {
  cat("First 3 sample entries:\n")
  for (i in 1:min(3, nrow(adaptive_data))) {
    cat("Position", adaptive_data$pos[i], "samples:", paste(adaptive_data$sample[[i]], collapse = ", "), "\n")
  }
  
  # Get all unique samples
  all_samples <- unique(unlist(adaptive_data$sample))
  cat("\nAll samples found:", paste(all_samples, collapse = ", "), "\n")
  cat("Number of unique samples:", length(all_samples), "\n")
}

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("This should help identify why adaptive_h4 results are all NA.\n")
