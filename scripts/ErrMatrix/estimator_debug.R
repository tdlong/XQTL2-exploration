#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Simple debug script: reproduce production results for one position
# Usage: Rscript scripts/ErrMatrix/estimator_debug.R <dump_rds>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/ErrMatrix/estimator_debug.R <dump_rds>")
}

dump_rds <- args[1]
payload <- readRDS(dump_rds)

# Extract
testing_position <- payload$testing_position
sample_name <- payload$sample_name
df_window <- payload$df_window
founders <- payload$founders
h_cutoff <- payload$h_cutoff
method <- payload$method
chr <- payload$chr

cat("Loaded payload:\n")
cat(" chr:", chr, " pos:", testing_position, " sample:", sample_name, "\n")
cat(" window SNPs:", nrow(df_window), " founders:", length(founders), "\n\n")

# Run the exact same logic as production (copy from est_haps_var)
cat("Running debug estimator logic...\n")

window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
accumulated_constraints <- NULL
accumulated_constraint_values <- NULL

V <- matrix(NA_real_, length(founders), length(founders))
rownames(V) <- founders; colnames(V) <- founders

for (window_size in window_sizes) {
  window_start <- max(1, testing_position - window_size/2)
  window_end <- testing_position + window_size/2
  window_data <- df_window %>% dplyr::filter(POS >= window_start & POS <= window_end)
  
  if (nrow(window_data) < 10) next
  if (!all(c(founders, sample_name) %in% names(window_data))) next

  founder_matrix <- window_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
  sample_freqs <- window_data %>% dplyr::pull(!!sample_name)
  complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  
  if (nrow(founder_matrix_clean) < 10) next

  founder_dist <- stats::dist(t(founder_matrix_clean))
  hclust_result <- stats::hclust(founder_dist, method = "complete")
  groups <- stats::cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))

  # Pooled covariance calculation (copy from production)
  gid <- sort(unique(groups))
  k <- length(gid)
  if (k == 0) next
  if (k == 1) {
    cov_pool <- matrix(1, 1, 1)
    pool_members <- list(which(groups == gid[1]))
  } else {
    A_pool <- matrix(0, nrow(founder_matrix_clean), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) founder_matrix_clean[, mem] else rowMeans(founder_matrix_clean[, mem, drop=FALSE])
    }
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(limSolve::lsei(A=A_pool, B=sample_freqs_clean, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar) && is.matrix(fit$covar) && nrow(fit$covar) == k && ncol(fit$covar) == k) {
      cov_pool <- fit$covar
      pool_members <- members
    } else {
      XtX <- crossprod(A_pool)
      xhat <- tryCatch(solve(XtX, crossprod(A_pool, sample_freqs_clean)), error=function(e) MASS::ginv(XtX) %*% crossprod(A_pool, sample_freqs_clean))
      r <- sample_freqs_clean - as.numeric(A_pool %*% xhat)
      sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
      cov_matrix <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
      if (!is.matrix(cov_matrix) || any(!is.finite(cov_matrix))) cov_matrix <- diag(k)
      cov_pool <- sigma2 * cov_matrix
      pool_members <- members
    }
  }

  # Update V matrix
  if (is.matrix(cov_pool) && nrow(cov_pool) == length(pool_members)) {
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        for (fi in mi) {
          for (fj in mj) {
            V[fi, fj] <- cov_pool[i_idx, j_idx]
          }
        }
      }
    }
  }

  # Check if all founders are separated
  if (n_groups == length(founders)) {
    cat("All founders separated at window", window_size, "bp\n")
    break
  }
}

# Build final result
groups_final <- unique(groups)
names_final <- founders[groups_final]
haps_final <- rep(1/length(groups_final), length(groups_final))

result <- list(
  Groups = groups_final,
  Haps = haps_final,
  Err = V,
  Names = names_final
)

cat("\n=== DEBUG RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse=","), "\n")
cat("Haplotypes:\n")
print(result$Haps)
cat("Error matrix:\n")
print(result$Err)
cat("Names:", paste(result$Names, collapse=","), "\n")

cat("\nDone.\n")


