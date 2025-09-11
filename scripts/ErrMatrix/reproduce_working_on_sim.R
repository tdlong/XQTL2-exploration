#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# --- BEGIN: Working estimation function copied from production (no edits) ---
estimate_haplotypes_list_format <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  method <- match.arg(method)
  if (verbose >= 2) {
    cat(sprintf("=== ADAPTIVE WINDOW: h_cutoff=%g, pos=%s, sample=%s ===\n", 
                h_cutoff, format(pos, big.mark=","), sample_name))
  }
  window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
  final_result <- NULL
  final_n_groups <- 0
  final_window_size <- window_sizes[1]
  final_wide_data <- NULL
  previous_n_groups <- 0
  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL
  for (window_idx in seq_along(window_sizes)) {
    window_size <- window_sizes[window_idx]
    window_start <- max(1, pos - window_size/2)
    window_end <- pos + window_size/2
    if (verbose >= 2) {
      cat(sprintf("  Window %d: %s", window_idx, 
                  ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp"))))
    }
    window_data <- df3 %>%
      dplyr::filter(.data$POS >= window_start & .data$POS <= window_end & .data$name %in% c(founders, sample_name))
    if (nrow(window_data) == 0) next
    wide_data <- window_data %>%
      dplyr::select(.data$POS, .data$name, .data$freq) %>%
      tidyr::pivot_wider(names_from = .data$name, values_from = .data$freq)
    if (!all(c(founders, sample_name) %in% names(wide_data)) || nrow(wide_data) < 10) next
    founder_matrix <- wide_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
    sample_freqs <- wide_data %>% dplyr::pull(!!sample_name)
    complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    if (nrow(founder_matrix_clean) < 10) next
    final_window_size <- window_size
    final_wide_data <- wide_data
    founder_dist <- stats::dist(t(founder_matrix_clean))
    hclust_result <- stats::hclust(founder_dist, method = "complete")
    groups <- stats::cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) cat(sprintf(" - %d SNPs, %d groups", nrow(founder_matrix_clean), n_groups))
    if (window_idx > 1 && n_groups <= previous_n_groups) { if (verbose >= 2) cat(" - ✗ No improvement\n"); next }
    previous_n_groups <- n_groups
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(rep(1, n_founders), nrow = 1)
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    ✓ Added %d accumulated constraints from previous windows\n", nrow(accumulated_constraints)))
    } else { if (verbose >= 2) cat("    • No accumulated constraints yet (first meaningful window)\n") }
    tryCatch({
      result <- limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean, 
                               E = E, F = F, 
                               G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE)
      if (result$IsError == 0) {
        current_constraints <- NULL
        current_constraint_values <- NULL
        if (verbose >= 2) cat(sprintf(" - ✓ LSEI success, %d groups", n_groups))
        unique_clusters <- unique(groups)
        for (cluster_id in unique_clusters) {
          cluster_founders <- which(groups == cluster_id)
          if (length(cluster_founders) > 1) {
            constraint_row <- rep(0, n_founders); constraint_row[cluster_founders] <- 1
            group_freq <- sum(result$X[cluster_founders])
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, group_freq)
          } else {
            founder_freq <- result$X[cluster_founders]
            constraint_row <- rep(0, n_founders); constraint_row[cluster_founders] <- 1
            current_constraints <- rbind(current_constraints, constraint_row)
            current_constraint_values <- c(current_constraint_values, founder_freq)
          }
        }
        if (!is.null(current_constraints)) {
          accumulated_constraints <- current_constraints
          accumulated_constraint_values <- current_constraint_values
          if (verbose >= 2) cat(sprintf(" - %d constraints", nrow(current_constraints)))
        } else {
          accumulated_constraints <- NULL; accumulated_constraint_values <- NULL
        }
        final_result <- result; final_n_groups <- n_groups
        if (verbose >= 2) {
          cat(sprintf("    Stored LSEI result with %d groups\n", n_groups))
          if ("covar" %in% names(result)) cat("    Covariance matrix available in this result (covar field)\n")
          else cat("    No covariance matrix in this result\n")
        }
        if (n_groups == length(founders)) { if (verbose >= 2) cat(" - ✅ SUCCESS: All founders distinguished\n"); break }
      }
    }, error = function(e) { if (verbose >= 2) cat(sprintf("LSEI error: %s\n", e$message)) })
  }
  if (!is.null(final_result)) {
    founder_frequencies <- final_result$X; names(founder_frequencies) <- founders
    estimate_OK <- if (final_n_groups == length(founders)) 1 else 0
  } else {
    founder_frequencies <- rep(NA, length(founders)); names(founder_frequencies) <- founders; estimate_OK <- NA
    if (verbose >= 1) cat("❌ LSEI failed or insufficient SNPs - no haplotype frequencies available\n")
  }
  if (!is.null(final_result)) {
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      error_matrix <- final_result$covar
      if (verbose >= 2) {
        cat("Using LSEI covariance matrix (covar field)\n")
        cat("\nError matrix (×100, rounded to 2 decimals):\n"); print(round(error_matrix * 100, 2)); cat("Condition number:", round(kappa(error_matrix), 2), "\n")
      }
    } else {
      error_matrix <- matrix(NA, length(founders), length(founders)); rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
      if (verbose >= 2) { cat("No covariance matrix in LSEI result - using NA matrix\n") }
    }
  } else {
    error_matrix <- matrix(NA, length(founders), length(founders)); rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
    if (verbose >= 2) { cat("No LSEI result - using NA matrix\n") }
  }
  return(list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders))
}
# --- END: Working estimation function copied ---

# --- Simulate founders and a sample, then run the working estimator ---
set.seed(123)
n_founders <- 8
founders <- paste0("F", 1:n_founders)
sample_name <- "S1"

# Simulate enough SNPs to cover largest window
n_snps <- 600000L
POS <- seq_len(n_snps)

# Base founder and related patterns
A <- matrix(0, nrow = n_snps, ncol = n_founders)
A[, 1] <- rbinom(n_snps, 1, 0.5)
block_size <- 150L
n_blocks <- n_snps %/% block_size
flip_in_blocks <- function(col, flips_per_block) {
  for (b in seq_len(n_blocks)) {
    s <- (b - 1L) * block_size + 1L; e <- min(b * block_size, n_snps)
    idx <- s:e
    if (length(idx) > 0 && flips_per_block > 0) {
      flip_idx <- sample(idx, size = min(flips_per_block, length(idx)))
      col[flip_idx] <- 1 - col[flip_idx]
    }
  }
  col
}
A[, 2] <- flip_in_blocks(A[, 1], 2)
A[, 3] <- flip_in_blocks(A[, 1], 10)
A[, 4] <- rbinom(n_snps, 1, 0.5)
A[, 5] <- flip_in_blocks(A[, 4], 4)
A[, 6] <- rbinom(n_snps, 1, 0.5)
A[, 7] <- rbinom(n_snps, 1, 0.5)
A[, 8] <- rbinom(n_snps, 1, 0.5)
colnames(A) <- founders

# Sample as mixture of founders (allele frequency)
set.seed(42)
w_true <- runif(n_founders); w_true <- w_true / sum(w_true)
sample_freq <- as.numeric(A %*% w_true) + rnorm(n_snps, 0, 0.02)
sample_freq[sample_freq < 0] <- 0; sample_freq[sample_freq > 1] <- 1

# Build df3-like table (POS, name, freq)
df_founders <- as_tibble(A) %>% mutate(POS = POS) %>%
  pivot_longer(cols = all_of(founders), names_to = "name", values_to = "freq")
df_sample <- tibble(POS = POS, name = sample_name, freq = sample_freq)
df3 <- bind_rows(df_founders, df_sample)

# Run the working estimator at mid position
pos_mid <- POS[round(length(POS)/2)]
h_cutoff <- 4
cat("=== Reproducing working estimator on simulated data ===\n")
cat("Position:", pos_mid, " h_cutoff:", h_cutoff, "\n")
res <- estimate_haplotypes_list_format(
  pos = pos_mid,
  sample_name = sample_name,
  df3 = df3,
  founders = founders,
  h_cutoff = h_cutoff,
  method = "adaptive",
  verbose = 2
)

cat("\n--- RESULT ---\n")
cat("Groups:", paste(res$Groups, collapse=","), "\n")
cat("Haps:\n"); print(round(res$Haps, 4))
if (is.matrix(res$Err)) {
  cat("\nErr condition number:\n"); print(kappa(res$Err))
} else { cat("\nNo error matrix available.\n") }


