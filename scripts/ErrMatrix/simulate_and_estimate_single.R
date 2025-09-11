#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Read-only import of production functions
source("scripts/production/complete_haplotype_workflow.R")

# 1) Simulate founders with controlled distinguishability (h_cutoff in Euclidean space)
set.seed(123)
n_founders <- 8
founders <- paste0("F", 1:n_founders)
sample_name <- "S1"

# Use enough SNPs to satisfy production window sizes (up to 500kb)
n_snps <- 600000L
POS <- seq_len(n_snps)

# Base founder
A <- matrix(0, nrow = n_snps, ncol = n_founders)
A[, 1] <- rbinom(n_snps, 1, 0.5)

# Create related founders by sparse flips per 150 positions to trigger 150/300/750/1500kb-like separations
block_size <- 150L
n_blocks <- n_snps %/% block_size

flip_in_blocks <- function(col, flips_per_block) {
  for (b in seq_len(n_blocks)) {
    start <- (b - 1L) * block_size + 1L
    end   <- min(b * block_size, n_snps)
    idx   <- start:end
    if (length(idx) > 0 && flips_per_block > 0) {
      flip_idx <- sample(idx, size = min(flips_per_block, length(idx)))
      col[flip_idx] <- 1 - col[flip_idx]
    }
  }
  col
}

# F2: 2 flips/150; F3: 10 flips/150; F4: random; F5: 4 flips/150; F6-8: random
A[, 2] <- flip_in_blocks(A[, 1], 2)
A[, 3] <- flip_in_blocks(A[, 1], 10)
A[, 4] <- rbinom(n_snps, 1, 0.5)
A[, 5] <- flip_in_blocks(A[, 4], 4)
A[, 6] <- rbinom(n_snps, 1, 0.5)
A[, 7] <- rbinom(n_snps, 1, 0.5)
A[, 8] <- rbinom(n_snps, 1, 0.5)

colnames(A) <- founders

# 2) Simulate sample as mixture of founders (allele frequency per SNP)
set.seed(42)
w_true <- runif(n_founders)
w_true <- w_true / sum(w_true)
sample_freq <- as.numeric(A %*% w_true) + rnorm(n_snps, 0, 0.02)
sample_freq[sample_freq < 0] <- 0
sample_freq[sample_freq > 1] <- 1

# 3) Build df3-compatible tibble: columns POS, name, freq for founders + sample
df_founders <- as_tibble(A) %>%
  mutate(POS = POS) %>%
  pivot_longer(cols = all_of(founders), names_to = "name", values_to = "freq")

df_sample <- tibble(POS = POS, name = sample_name, freq = sample_freq)

df3 <- bind_rows(df_founders, df_sample)

# 4) Choose a mid-genome position so windows have data on both sides
pos_mid <- POS[round(length(POS)/2)]
h_cutoff <- 4

cat("=== Simulated estimation at position:", pos_mid, "(h_cutoff =", h_cutoff, ")\n")

res <- estimate_haplotypes_list_format(
  pos = pos_mid,
  sample_name = sample_name,
  df3 = df3,
  founders = founders,
  h_cutoff = h_cutoff,
  method = "adaptive",
  window_size_bp = NULL,
  chr = "chr2R",
  verbose = 2
)

cat("\n--- RESULT ---\n")
cat("Groups:", paste(res$Groups, collapse=","), "\n")
cat("Haps (first 8):\n")
print(round(res$Haps, 4))

if (is.matrix(res$Err)) {
  cat("\nError matrix condition number:\n")
  print(kappa(res$Err))
} else {
  cat("\nNo error matrix available.\n")
}


