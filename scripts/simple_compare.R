#!/usr/bin/env Rscript

library(dplyr)

# Load the two datasets with correct column names
new_estimates <- readRDS("process/JUICE/haplotype_results/snp_imputation_adaptive_h8_chr2R_DEBUG.RDS") %>% select(POS, imputed) %>% rename(new_estimate = imputed)
old_estimates <- readRDS("process/JUICE/haplotype_results/snp_imputation_adaptive_h8_chr2R.RDS") %>% filter(sample == "AJ_1_1") %>% select(pos, imputed) %>% rename(POS = pos, old_estimate = imputed)

# Left join and compare
comparison <- new_estimates %>% left_join(old_estimates, by = "POS") %>% mutate(diff = abs(new_estimate - old_estimate))

print(comparison)
cat("\nMean difference:", mean(comparison$diff, na.rm = TRUE), "\n")
cat("Max difference:", max(comparison$diff, na.rm = TRUE), "\n")
