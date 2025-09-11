#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Read-only: source production workflow to access the WORKING function
source("scripts/production/complete_haplotype_workflow.R")

# Keep a handle to the working function (DO NOT EDIT PRODUCTION)
estimate_haplotypes_list_format_prod <- estimate_haplotypes_list_format

# Editable copy with the same signature/return contract.
# Behavior:
# - If pos == -99: use prefix-based windows (1..150, 1..300, 1..750, 1..1500, 1..3000)
#   on the provided df3 (POS assumed 1..N) and run adaptive LSEI exactly like production.
# - Else: delegate to production function unchanged.
estimate_haplotypes_list_format_sim <- function(pos, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {

  # Delegate to production when not in simulated mode
  if (!isTRUE(all.equal(pos, -99))) {
    return(estimate_haplotypes_list_format_prod(pos, sample_name, df3, founders, h_cutoff,
                                                method, window_size_bp, chr, verbose))
  }

  if (verbose >= 2) {
    cat(sprintf("=== SIM MODE (-99): h_cutoff=%g, sample=%s ===\n", h_cutoff, sample_name))
  }

  window_sizes <- c(150, 300, 750, 1500, 3000)
  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founders))

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  # Ensure POS-ordered tidy input (POS, name, freq)
  df3 <- df3 %>% dplyr::arrange(.data$POS)

  for (window_size in window_sizes) {
    window_data <- df3 %>%
      dplyr::filter(.data$POS >= 1 & .data$POS <= window_size & .data$name %in% c(founders, sample_name))
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

    # Clustering (Euclidean, complete linkage), cutree at h_cutoff
    founder_dist <- stats::dist(t(founder_matrix_clean))
    hclust_result <- stats::hclust(founder_dist, method = "complete")
    groups <- stats::cutree(hclust_result, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) cat(sprintf("  Prefix %d SNPs → %d groups\n", nrow(founder_matrix_clean), n_groups))

    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("    No improvement; continue\n")
      next
    }
    previous_n_groups <- n_groups

    # Constraints accumulate like production
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    Added %d accumulated constraints\n", nrow(accumulated_constraints)))
    }

    res <- tryCatch(
      limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean,
                     E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("    LSEI failed; continue\n")
      next
    }

    # Build/accumulate group constraints (pool groups, lock singles)
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

  # Return EXACT same structure as production
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
    error_matrix <- matrix(NA, length(founders), length(founders))
    rownames(error_matrix) <- founders; colnames(error_matrix) <- founders
  }

  list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders)
}

# Optional demo when run directly
if (sys.nframe() == 0) {
  # No demo here. Use your external wrapper to simulate and call this file's API.
}

# Simple alias wrapper for clarity in callers
estimate_haplotypes_dev <- function(pos, sample_name, df3, founders, h_cutoff,
                                   method = "adaptive", window_size_bp = NULL,
                                   chr = "chr2R", verbose = 0) {
  estimate_haplotypes_list_format_sim(pos, sample_name, df3, founders, h_cutoff,
                                      method, window_size_bp, chr, verbose)
}

# Wrapper function: simulate founders + sample → df3, then call the modified estimator
run_simulation_wrapper <- function(n_snps = 3000L, h_cutoff = 4, verbose = 1) {
  set.seed(123)
  n_founders <- 8
  founders <- paste0("F", seq_len(n_founders))
  sample_name <- "S1"
  POS <- seq_len(n_snps)

  # Simulate founders (controlled distinguishability)
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
  samp <- as.numeric(A %*% w) + rnorm(n_snps, 0, 0.02)
  samp[samp < 0] <- 0; samp[samp > 1] <- 1

  # Build df3-like tibble
  df_founders <- tibble::as_tibble(A) %>% dplyr::mutate(POS = POS) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(founders), names_to = "name", values_to = "freq")
  df_sample <- tibble::tibble(POS = POS, name = sample_name, freq = samp)
  df3 <- dplyr::bind_rows(df_founders, df_sample)

  # Call the modified estimator in prefix mode (pos = -99)
  res <- estimate_haplotypes_list_format_sim(
    pos = -99, sample_name = sample_name, df3 = df3, founders = founders,
    h_cutoff = h_cutoff, method = "adaptive", verbose = verbose
  )

  if (verbose >= 1) {
    cat("Groups:", paste(res$Groups, collapse=","), "\n")
    print(round(res$Haps, 4))
    if (is.matrix(res$Err)) cat("Err kappa:", kappa(res$Err), "\n")
  }

  invisible(res)
}


