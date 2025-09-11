#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# 1) Source production workflow (read-only). We'll shadow just the one function locally.
source("scripts/production/complete_haplotype_workflow.R")

# Stash original function so we can call it when pos != -99
original_estimate_haplotypes_list_format <- estimate_haplotypes_list_format

# 2) Define a wrapper that preserves the exact signature/return, but supports pos == -99
# If pos == -99, we interpret df3 as a simulated 1..N POS and grow prefix windows
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {

  # If not the special simulated mode, delegate to production function unchanged
  if (!is.numeric(pos) || length(pos) != 1 || (!is.na(pos) && pos > 0)) {
    return(original_estimate_haplotypes_list_format(pos, sample_name, df3, founders, h_cutoff,
                                                    method, window_size_bp, chr, verbose))
  }
  if (!identical(pos, -99)) {
    return(original_estimate_haplotypes_list_format(pos, sample_name, df3, founders, h_cutoff,
                                                    method, window_size_bp, chr, verbose))
  }

  # Simulated prefix-window mode (pos == -99): grow windows 1..150, 1..300, ... over POS 1..N
  if (verbose >= 2) {
    cat(sprintf("=== PREFIX MODE (simulated): h_cutoff=%g, sample=%s ===\n", h_cutoff, sample_name))
  }

  window_sizes <- c(150, 300, 750, 1500, 3000)
  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founders))

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  # df3 should contain POS, name, freq for all founders and the sample
  # Ensure ordering by POS
  df3 <- df3 %>% arrange(.data$POS)

  for (window_size in window_sizes) {
    # Slice prefix
    window_data <- df3 %>%
      filter(.data$POS >= 1 & .data$POS <= window_size & .data$name %in% c(founders, sample_name))
    if (nrow(window_data) == 0) next

    wide_data <- window_data %>%
      select(.data$POS, .data$name, .data$freq) %>%
      pivot_wider(names_from = .data$name, values_from = .data$freq)

    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next

    founder_matrix <- wide_data %>% select(all_of(founders)) %>% as.matrix()
    sample_freqs <- wide_data %>% pull(!!sample_name)

    complete_rows <- complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    if (nrow(founder_matrix_clean) < 10) next

    # Clustering
    founder_dist <- dist(t(founder_matrix_clean))
    hclust_result <- hclust(founder_dist, method = "complete")
    groups <- cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) cat(sprintf("  Prefix %d SNPs → %d groups\n", nrow(founder_matrix_clean), n_groups))

    # Require improvement over previous step
    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("    No improvement; continue\n")
      next
    }
    previous_n_groups <- n_groups

    # Constraints
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    Added %d accumulated constraints\n", nrow(accumulated_constraints)))
    }

    # Solve LSEI
    res <- tryCatch(
      limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean,
                     E = E, F = F,
                     G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("    LSEI failed; continue\n")
      next
    }

    # Accumulate group constraints like production
    current_constraints <- NULL
    current_constraint_values <- NULL
    for (cid in unique(groups)) {
      idx <- which(groups == cid)
      row <- rep(0, n_founders); row[idx] <- 1
      current_constraints <- rbind(current_constraints, row)
      current_constraint_values <- c(current_constraint_values, sum(res$X[idx]))
    }
    if (!is.null(current_constraints)) {
      accumulated_constraints <- current_constraints
      accumulated_constraint_values <- current_constraint_values
    } else {
      accumulated_constraints <- NULL
      accumulated_constraint_values <- NULL
    }

    final_result <- res
    final_n_groups <- n_groups
    if (verbose >= 2) cat("    Stored result\n")
    if (n_groups == length(founders)) { if (verbose >= 2) cat("    All founders distinguished; stop\n"); break }
  }

  # Build outputs identical to production
  if (!is.null(final_result)) {
    founder_frequencies <- final_result$X; names(founder_frequencies) <- founders
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      error_matrix <- final_result$covar
    } else {
      error_matrix <- matrix(NA, length(founders), length(founders))
      rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
    }
  } else {
    founder_frequencies <- rep(NA_real_, length(founders)); names(founder_frequencies) <- founders
    error_matrix <- matrix(NA, length(founders), length(founders)); rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
  }

  return(list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders))
}

# 3) Simulate founders and sample, make df3 (POS 1..3000), then call with pos = -99
set.seed(123)
n_founders <- 8
founders <- paste0("F", seq_len(n_founders))
sample_name <- "S1"
n_snps <- 3000L
POS <- seq_len(n_snps)

# Simulate founders with controlled distinguishability
A <- matrix(0, nrow = n_snps, ncol = n_founders)
A[, 1] <- rbinom(n_snps, 1, 0.5)
block <- 150L; n_blocks <- n_snps %/% block
flip_blocks <- function(col, k) {
  for (b in seq_len(n_blocks)) {
    s <- (b - 1L) * block + 1L; e <- min(b * block, n_snps)
    idx <- s:e
    if (k > 0) { flip <- sample(idx, size = min(k, length(idx))); col[flip] <- 1 - col[flip] }
  }
  col
}
A[, 2] <- flip_blocks(A[, 1], 2)
A[, 3] <- flip_blocks(A[, 1], 10)
A[, 4] <- rbinom(n_snps, 1, 0.5)
A[, 5] <- flip_blocks(A[, 4], 4)
A[, 6] <- rbinom(n_snps, 1, 0.5)
A[, 7] <- rbinom(n_snps, 1, 0.5)
A[, 8] <- rbinom(n_snps, 1, 0.5)
colnames(A) <- founders

set.seed(42)
w <- runif(n_founders); w <- w / sum(w)
sample_freq <- as.numeric(A %*% w) + rnorm(n_snps, 0, 0.02)
sample_freq[sample_freq < 0] <- 0; sample_freq[sample_freq > 1] <- 1

# Build df3-like tibble: POS, name, freq for founders and sample
df_founders <- as_tibble(A) %>% mutate(POS = POS) %>% pivot_longer(cols = all_of(founders), names_to = "name", values_to = "freq")
df_sample <- tibble(POS = POS, name = sample_name, freq = sample_freq)
df3 <- bind_rows(df_founders, df_sample)

cat("=== Running wrapper in simulated prefix mode (pos = -99) ===\n")
res <- estimate_haplotypes_list_format(
  pos = -99, sample_name = sample_name, df3 = df3, founders = founders, h_cutoff = 4,
  method = "adaptive", window_size_bp = NULL, chr = "chr2R", verbose = 2
)

cat("\n--- RESULT ---\n")
cat("Groups:", paste(res$Groups, collapse=","), "\n")
cat("Haps:\n"); print(round(res$Haps, 4))
if (is.matrix(res$Err)) { cat("\nErr condition number:\n"); print(kappa(res$Err)) } else { cat("\nNo error matrix available.\n") }


