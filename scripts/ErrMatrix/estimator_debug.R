#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
})

# Standalone estimator runner with diagnostics
# Usage:
#   Rscript scripts/ErrMatrix/estimator_debug.R <dump_rds>
# The RDS must be produced by BASE_VAR_WIDE.R with --dump-estimator-input

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

# Minimal copy of estimator's core up to pooled_cov computation with diagnostics

window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)

diagnostics <- list()

accumulated_constraints <- NULL
accumulated_constraint_values <- NULL

V <- matrix(NA_real_, length(founders), length(founders))
rownames(V) <- founders; colnames(V) <- founders

pooled_cov <- function(A_full, y_full, groups_vec){
  gid <- sort(unique(groups_vec))
  k <- length(gid)
  if (k == 0) return(list(cov=matrix(NA, 0, 0), members=list(), w=numeric(0)))
  if (k == 1) return(list(cov=matrix(1, 1, 1), members=list(1), w=1))
  A_pool <- matrix(0, nrow(A_full), k)
  members <- vector("list", k)
  for (ii in seq_along(gid)){
    mem <- which(groups_vec == gid[ii])
    members[[ii]] <- mem
    A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
  }
  E <- matrix(1, 1, k); F <- 1
  fit <- tryCatch(limSolve::lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
  if (!is.null(fit) && !is.null(fit$covar) && is.matrix(fit$covar) && nrow(fit$covar) == k && ncol(fit$covar) == k) {
    return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
  }
  XtX <- crossprod(A_pool)
  xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) MASS::ginv(XtX) %*% crossprod(A_pool, y_full))
  r <- y_full - as.numeric(A_pool %*% xhat)
  sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
  cov_matrix <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
  if (!is.matrix(cov_matrix) || any(!is.finite(cov_matrix))) cov_matrix <- diag(k)
  list(cov = sigma2 * cov_matrix, members=members, w=as.numeric(xhat))
}

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

  pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
  cov_pool <- pc$cov; pool_members <- pc$members

  diag_sum <- sum(abs(diag(cov_pool)), na.rm = TRUE)
  off <- cov_pool; diag(off) <- 0
  off_sum <- sum(abs(off), na.rm = TRUE)
  kappa_val <- tryCatch(kappa(cov_pool), error=function(e) NA_real_)
  eig_vals <- tryCatch(eigen(cov_pool, only.values=TRUE)$values, error=function(e) NA_real_)
  eig_min <- if (all(is.na(eig_vals))) NA_real_ else min(Re(eig_vals), na.rm=TRUE)
  eig_max <- if (all(is.na(eig_vals))) NA_real_ else max(Re(eig_vals), na.rm=TRUE)

  diagnostics[[as.character(window_size)]] <- list(
    window_size = window_size,
    n_snps = nrow(founder_matrix_clean),
    n_groups = n_groups,
    pool_members_sizes = purrr::map_int(pool_members, length),
    cov_dim = dim(cov_pool),
    kappa = kappa_val,
    eig_min = eig_min,
    eig_max = eig_max,
    diag_sum = diag_sum,
    off_sum = off_sum
  )
}

cat("Diagnostics by window size (bp):\n")
purrr::iwalk(diagnostics, function(d, nm){
  cat(sprintf("%6s: n_snps=%-5d n_groups=%-2d cov=%sx%s kappa=%.3e eig=[%.3e, %.3e] diag=%.3e off=%.3e\n",
              nm, d$n_snps, d$n_groups, d$cov_dim[1], d$cov_dim[2], d$kappa, d$eig_min, d$eig_max, d$diag_sum, d$off_sum))
})

cat("\nDone.\n")


