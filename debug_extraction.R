#!/usr/bin/env Rscript

library(tidyverse)

# Load files
orig <- readRDS('haplotype_results_list_format/adaptive_window_h4_results_chr3R.RDS')
resh <- readRDS('haplotype_results_list_format/adapt_h4/R.haps.chr3R.out.rds')

# Unnest reshaped
dash_u <- resh %>%
  tidyr::unnest(c(sample, Groups, Haps, Err, Names))

# Get first row for comparison
orig_row <- orig[1,]
dash_row <- dash_u[1,]

cat("=== ORIGINAL ROW ===\n")
cat("Names class:", class(orig_row$Names), "length:", length(orig_row$Names), "\n")
cat("Names content:", paste(orig_row$Names, collapse=","), "\n")
cat("Err class:", class(orig_row$Err), "is.matrix:", is.matrix(orig_row$Err), "\n")
cat("Err dim:", paste(dim(orig_row$Err), collapse="x"), "\n")
cat("Err rownames:", paste(rownames(orig_row$Err), collapse=","), "\n")

cat("\n=== RESHAPED ROW (after unnesting) ===\n")
cat("Names class:", class(dash_row$Names), "length:", length(dash_row$Names), "\n")
cat("Names content:", paste(dash_row$Names, collapse=","), "\n")
cat("Err class:", class(dash_row$Err), "is.matrix:", is.matrix(dash_row$Err), "\n")
cat("Err dim:", paste(dim(dash_row$Err), collapse="x"), "\n")
cat("Err rownames:", paste(rownames(dash_row$Err), collapse=","), "\n")

# Test the join
key_cols <- c("CHROM","pos","sample")
orig_key <- orig %>% select(all_of(key_cols)) %>% mutate(in_orig = TRUE)
dash_key <- dash_u %>% select(all_of(key_cols)) %>% mutate(in_resh = TRUE)
keys <- full_join(orig_key, dash_key, by = key_cols)
common <- keys %>% filter(in_orig == TRUE, in_resh == TRUE) %>% select(all_of(key_cols))

# Join payloads
dfc <- common %>%
  left_join(orig, by = key_cols) %>%
  rename(Groups_o = Groups, Haps_o = Haps, Err_o = Err, Names_o = Names) %>%
  left_join(dash_u, by = key_cols) %>%
  rename(Groups_r = Groups, Haps_r = Haps, Err_r = Err, Names_r = Names)

cat("\n=== JOINED ROW ===\n")
row <- dfc[1,]
cat("Names_o class:", class(row$Names_o), "length:", length(row$Names_o), "\n")
cat("Names_r class:", class(row$Names_r), "length:", length(row$Names_r), "\n")
cat("Err_o class:", class(row$Err_o), "is.matrix:", is.matrix(row$Err_o), "\n")
cat("Err_r class:", class(row$Err_r), "is.matrix:", is.matrix(row$Err_r), "\n")
cat("Err_r rownames length:", length(rownames(row$Err_r)), "\n")
