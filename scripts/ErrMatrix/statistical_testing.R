#!/usr/bin/env Rscript

# Statistical testing for haplotype estimation methods
# Usage: Rscript statistical_testing.R <data_dir> <chromosome> <start_pos> <end_pos> <design_file>

library(tidyverse)
library(limSolve)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript statistical_testing.R <data_dir> <chromosome> <start_pos> <end_pos> <design_file>")
}

mydir <- args[1]
chromosome <- args[2]
start_pos <- as.numeric(args[3])
end_pos <- as.numeric(args[4])
design_file <- args[5]

cat("=== STATISTICAL TESTING PARAMETERS ===\n")
cat("Data directory:", mydir, "\n")
cat("Chromosome:", chromosome, "\n")
cat("Position range:", start_pos, "-", end_pos, "\n")
cat("Design file:", design_file, "\n")

# Read design file
design.df = read.table(design_file, header=T, stringsAsFactors=F)
cat("Design matrix:\n")
print(design.df)

# Read haplotype data
input_file <- file.path(mydir, paste0("R.haps.", chromosome, ".out.rds"))
cat("Input file:", input_file, "\n")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

xx1 = readRDS(input_file)

# Wald test function
wald.test3 = function(p1, p2, covar1, covar2) {
  # p1, p2: haplotype frequencies for treatments 1 and 2
  # covar1, covar2: covariance matrices for treatments 1 and 2
  
  # Calculate difference
  diff = p1 - p2
  
  # Calculate combined covariance
  covar_combined = covar1 + covar2
  
  # Calculate Wald statistic
  tryCatch({
    # Use generalized inverse for numerical stability
    covar_inv = ginv(covar_combined)
    wald_stat = t(diff) %*% covar_inv %*% diff
    p_value = 1 - pchisq(wald_stat, df = length(diff))
    log10p = -log10(p_value)
    return(log10p)
  }, error = function(e) {
    return(NA)
  })
}

# Main analysis function
doscan2 = function(data, CHROM) {
  # Extract sample names and create design matrix
  sample_names <- data$sample[[1]]
  design_df <- data.frame(
    sample = sample_names,
    bam = sample_names,
    TRT = ifelse(str_detect(sample_names, "_W_"), "W", "Z"),
    REP = str_extract(sample_names, "Rep\\d+"),
    stringsAsFactors = FALSE
  ) %>%
    left_join(design.df, by = "bam") %>%
    filter(!is.na(TRT))
  
  # Filter for male samples only
  design_df <- design_df %>% filter(str_detect(sample, "_M"))
  
  if (nrow(design_df) == 0) {
    return(list(Wald_log10p = NA))
  }
  
  # Separate by treatment
  C_samples <- design_df %>% filter(TRT == "C") %>% pull(sample)
  Z_samples <- design_df %>% filter(TRT == "Z") %>% pull(sample)
  
  if (length(C_samples) == 0 || length(Z_samples) == 0) {
    return(list(Wald_log10p = NA))
  }
  
  # Extract haplotype frequencies and error matrices
  hap_freqs <- data$Haps[[1]]
  err_matrices <- data$Err[[1]]
  
  # Get indices for each treatment
  C_indices <- which(sample_names %in% C_samples)
  Z_indices <- which(sample_names %in% Z_samples)
  
  # Calculate average haplotype frequencies per treatment
  p1 <- colMeans(hap_freqs[C_indices, , drop = FALSE])  # Treatment C
  p2 <- colMeans(hap_freqs[Z_indices, , drop = FALSE])  # Treatment Z
  
  # Calculate average error matrices per treatment
  covar1 <- Reduce(`+`, err_matrices[C_indices]) / length(C_indices)  # Treatment C
  covar2 <- Reduce(`+`, err_matrices[Z_indices]) / length(Z_indices)  # Treatment Z
  
  # Perform Wald test
  Wald_log10p <- wald.test3(p1, p2, covar1, covar2)
  
  # Calculate treatment differences and error variances
  hap_diff <- p1 - p2
  err_var_C <- diag(covar1)
  err_var_Z <- diag(covar2)
  err_diff <- err_var_C - err_var_Z
  
  return(list(
    Wald_log10p = Wald_log10p,
    hap_freqs_C = list(p1),
    hap_freqs_Z = list(p2), 
    hap_diff = list(hap_diff),
    err_var_C = list(err_var_C),
    err_var_Z = list(err_var_Z),
    err_diff = list(err_diff)
  ))
}

# Run analysis
Nfounders=length(xx1$Groups[[1]][[1]])

bb1 = xx1 %>%
  filter(pos >= start_pos & pos <= end_pos) %>%
  group_by(CHROM,pos) %>%
  nest() %>%
  mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
  unnest_wider(out)

bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)

# Display Wald test results
cat("\n=== WALD TEST RESULTS ===\n")
print(bb2 %>% select(pos, Wald_log10p), n = Inf)

# HAPLOTYPE CHANGES (NUMERATOR)
cat("\n=== HAPLOTYPE CHANGES BETWEEN ADJACENT POSITIONS ===\n")

hap_analysis <- bb2 %>%
  filter(!is.na(Wald_log10p)) %>%
  arrange(pos) %>%
  mutate(hap_diff = map(hap_diff, ~ .x[[1]])) %>%
  select(pos, hap_diff)

hap_changes <- hap_analysis %>%
  mutate(
    hap_diff_prev = lag(hap_diff),
    hap_change = map2(hap_diff, hap_diff_prev, ~ {
      if (is.null(.y)) return(NULL)
      abs(.x - .y)
    })
  ) %>%
  filter(!is.null(hap_diff_prev)) %>%
  select(pos, hap_change)

if (nrow(hap_changes) > 0) {
  cat("Changes in haplotype differences (×1000):\n")
  cat("Position     B1    B2    B3    B4    B5    B6    B7   AB8\n")
  for(i in 1:nrow(hap_changes)) {
    if (!is.null(hap_changes$hap_change[[i]])) {
      changes <- round(hap_changes$hap_change[[i]] * 1000, 0)
      cat(sprintf("%9.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n", 
                  hap_changes$pos[i], changes[1], changes[2], changes[3], changes[4], 
                  changes[5], changes[6], changes[7], changes[8]))
    }
  }
  
  all_changes <- unlist(hap_changes$hap_change)
  cat(sprintf("\nHaplotype change summary (×1000): Mean=%.1f, SD=%.1f, Max=%.1f\n", 
              mean(all_changes) * 1000, sd(all_changes) * 1000, max(all_changes) * 1000))
} else {
  cat("No haplotype changes calculated\n")
}

# ERROR VARIANCE CHANGES (DENOMINATOR)
cat("\n=== ERROR VARIANCE CHANGES BETWEEN ADJACENT POSITIONS ===\n")

err_analysis <- bb2 %>%
  filter(!is.na(Wald_log10p)) %>%
  arrange(pos) %>%
  mutate(
    err_var_C = map(err_var_C, ~ .x[[1]]),
    err_var_Z = map(err_var_Z, ~ .x[[1]]),
    sqrt_avg_err_var = map2(err_var_C, err_var_Z, ~ sqrt((.x + .y) / 2))
  ) %>%
  select(pos, sqrt_avg_err_var)

err_changes <- err_analysis %>%
  mutate(
    err_prev = lag(sqrt_avg_err_var),
    err_change = map2(sqrt_avg_err_var, err_prev, ~ {
      if (is.null(.y)) return(NULL)
      abs(.x - .y)
    })
  ) %>%
  filter(!is.null(err_prev)) %>%
  select(pos, err_change)

if (nrow(err_changes) > 0) {
  cat("Changes in sqrt(average variance) (×1000):\n")
  cat("Position     B1    B2    B3    B4    B5    B6    B7   AB8\n")
  for(i in 1:nrow(err_changes)) {
    if (!is.null(err_changes$err_change[[i]])) {
      changes <- round(err_changes$err_change[[i]] * 1000, 0)
      cat(sprintf("%9.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n", 
                  err_changes$pos[i], changes[1], changes[2], changes[3], changes[4], 
                  changes[5], changes[6], changes[7], changes[8]))
    }
  }
  
  all_err_changes <- unlist(err_changes$err_change)
  cat(sprintf("\nError variance change summary (×1000): Mean=%.1f, SD=%.1f, Max=%.1f\n", 
              mean(all_err_changes) * 1000, sd(all_err_changes) * 1000, max(all_err_changes) * 1000))
} else {
  cat("No error variance changes calculated\n")
}