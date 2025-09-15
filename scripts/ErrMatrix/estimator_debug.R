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
  stop("Usage: Rscript scripts/ErrMatrix/estimator_debug.R <dump_rds> [--v1] [--v2]")
}

dump_rds <- args[1]
verbosity <- 0
if ("--v1" %in% args) verbosity <- 1
if ("--v2" %in% args) verbosity <- 2
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

fmt_cell_signed <- function(x, diag_cell){
  if (is.na(x)) return("  NA ")
  if (x == 0) {
    signc <- if (diag_cell) " " else "+"
    return(paste0(signc, "  0 "))
  }
  signc <- if (diag_cell) " " else if (x < 0) "-" else "+"
  ax <- abs(x)
  expo <- floor(log10(ax))
  mant <- ax / (10^expo)
  m2 <- as.integer(round(mant * 10))
  if (expo < -9) expo <- -9
  if (expo > 9) expo <- 9
  paste0(signc, sprintf("%02d%+1d", m2, as.integer(expo)))
}

print_V_compact <- function(Vmat){
  cat("V (signed covariances, sci 2d+exp; diag no sign):\n")
  for (i in seq_len(nrow(Vmat))){
    row <- character(ncol(Vmat))
    for (j in seq_len(ncol(Vmat))){
      row[j] <- fmt_cell_signed(Vmat[i, j], diag_cell = (i == j))
    }
    cat(paste(row, collapse=" "), "\n", sep="")
  }
}

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
  n_snps_raw <- nrow(window_data)
  if (n_snps_raw < 10) next
  if (!all(c(founders, sample_name) %in% names(window_data))) next

  founder_matrix <- window_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
  sample_freqs <- window_data %>% dplyr::pull(!!sample_name)
  complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
  founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
  sample_freqs_clean <- sample_freqs[complete_rows]
  n_snps_clean <- nrow(founder_matrix_clean)
  if (n_snps_clean < 10) next

  founder_dist <- stats::dist(t(founder_matrix_clean))
  hclust_result <- stats::hclust(founder_dist, method = "complete")
  groups <- stats::cutree(hclust_result, h = h_cutoff)
  n_groups <- length(unique(groups))

  pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
  cov_pool <- pc$cov; pool_members <- pc$members

  # Render groups like production (e.g., [12]3[45]678)
  grp_ids <- sort(unique(groups))
  group_strings <- {
    parts <- character(0)
    for (gid in grp_ids) {
      idx <- which(groups == gid)
      if (length(idx) > 1) {
        parts <- c(parts, paste0("[", paste(idx, collapse=""), "]"))
      } else {
        parts <- c(parts, paste0(idx))
      }
    }
    paste(parts, collapse = "")
  }

  # Build a one-step V snapshot for singleton-vs-singleton covariances
  V_step <- matrix(NA_real_, length(founders), length(founders))
  rownames(V_step) <- founders; colnames(V_step) <- founders
  if (is.matrix(cov_pool)) {
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        if (length(mi)==1L && length(mj)==1L){
          fi <- founders[mi]; fj <- founders[mj]
          V_step[fi, fj] <- cov_pool[i_idx, j_idx]
          V_step[fj, fi] <- cov_pool[j_idx, i_idx]
        }
      }
    }
  }

  if (verbosity >= 1) {
    cat(sprintf("\n[%6s bp] groups: %s  (n_snps=%d)\n", format(window_size, trim=TRUE), group_strings, n_snps_clean))
    if (is.matrix(cov_pool)) {
      cat("V snapshot (singletons only):\n")
      print_V_compact(V_step)
    } else {
      cat("V snapshot: cov_pool unavailable (non-matrix)\n")
    }
  }

  if (verbosity >= 2) {
    # Detailed prints
    cat("pool_members (by founder indices): ", paste(vapply(pool_members, function(x) paste(x, collapse=","), ""), collapse=" | "), "\n", sep="")
    # Show SVD singular values of centered A_pool
    if (exists("A_pool")) {
      svd_vals <- tryCatch(svd(scale(A_pool, center=TRUE, scale=FALSE))$d, error=function(e) NA_real_)
      cat("A_pool singular values:", paste(signif(svd_vals, 4), collapse=", "), "\n")
    }
    # Top correlated founder pairs
    if (exists("founder_cor") && !is.null(founder_cor)) {
      off_mask <- row(founder_cor) < col(founder_cor)
      if (any(off_mask, na.rm=TRUE)) {
        vals <- abs(founder_cor[off_mask])
        ord <- order(vals, decreasing=TRUE, na.last=NA)
        top_n <- min(5, length(ord))
        if (top_n > 0) {
          cat("Top founder |cor| pairs:\n")
          for (ii in seq_len(top_n)) {
            # find the i,j for ord[ii]
            idx <- which(off_mask, arr.ind=TRUE)[ord[ii], ]
            i <- idx[1]; j <- idx[2]
            cat(sprintf("  %s-%s: %.4f\n", colnames(founder_cor)[i], colnames(founder_cor)[j], abs(founder_cor[i,j])))
          }
        }
      }
    }
    # cov_pool matrix print (full)
    if (is.matrix(cov_pool)) {
      cat("cov_pool (full):\n")
      print(signif(cov_pool, 4))
    }
  }

  # Additional diagnostics
  # 1) NA filtering stats and ranges
  na_dropped <- n_snps_raw - n_snps_clean
  sample_range <- c(min(sample_freqs_clean, na.rm=TRUE), max(sample_freqs_clean, na.rm=TRUE))
  founder_ranges <- apply(founder_matrix_clean, 2, function(col) c(min(col, na.rm=TRUE), max(col, na.rm=TRUE)))
  # 2) Pairwise correlations among founders (on clean SNPs)
  founder_cor <- tryCatch(cor(founder_matrix_clean, use="pairwise.complete.obs"), error=function(e) NULL)
  founder_cor_offdiag_max <- if (!is.null(founder_cor)) max(abs(founder_cor[row(founder_cor)!=col(founder_cor)]), na.rm=TRUE) else NA_real_
  founder_cor_offdiag_min <- if (!is.null(founder_cor)) min(founder_cor[row(founder_cor)!=col(founder_cor)], na.rm=TRUE) else NA_real_
  # 3) Build pooled design and analyze SVD/Rank/Correlation
  gid <- sort(unique(groups))
  k <- length(gid)
  A_pool <- matrix(0, nrow(founder_matrix_clean), k)
  for (ii in seq_along(gid)){
    mem <- which(groups == gid[ii])
    A_pool[, ii] <- if (length(mem) == 1L) founder_matrix_clean[, mem] else rowMeans(as.matrix(founder_matrix_clean[, mem, drop=FALSE]))
  }
  sv <- tryCatch(svd(scale(A_pool, center=TRUE, scale=FALSE))$d, error=function(e) NA_real_)
  rank_est <- if (all(is.na(sv))) NA_integer_ else sum(sv > max(1e-12, 1e-8*max(sv, na.rm=TRUE)), na.rm=TRUE)
  ap_cor <- tryCatch(cor(A_pool), error=function(e) NULL)
  ap_cor_offdiag_max <- if (!is.null(ap_cor)) max(abs(ap_cor[row(ap_cor)!=col(ap_cor)]), na.rm=TRUE) else NA_real_
  # 4) XtX eig/cond and residual variance
  XtX <- crossprod(A_pool)
  xhat <- tryCatch(solve(XtX, crossprod(A_pool, sample_freqs_clean)), error=function(e) MASS::ginv(XtX) %*% crossprod(A_pool, sample_freqs_clean))
  r <- sample_freqs_clean - as.numeric(A_pool %*% xhat)
  sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
  XtX_eigs <- tryCatch(eigen(XtX, only.values=TRUE)$values, error=function(e) NA_real_)
  XtX_eig_min <- if (all(is.na(XtX_eigs))) NA_real_ else min(Re(XtX_eigs), na.rm=TRUE)
  XtX_eig_max <- if (all(is.na(XtX_eigs))) NA_real_ else max(Re(XtX_eigs), na.rm=TRUE)
  XtX_kappa <- tryCatch(kappa(XtX), error=function(e) NA_real_)

  diag_sum <- sum(abs(diag(cov_pool)), na.rm = TRUE)
  off <- cov_pool; diag(off) <- 0
  off_sum <- sum(abs(off), na.rm = TRUE)
  kappa_val <- tryCatch(kappa(cov_pool), error=function(e) NA_real_)
  eig_vals <- tryCatch(eigen(cov_pool, only.values=TRUE)$values, error=function(e) NA_real_)
  eig_min <- if (all(is.na(eig_vals))) NA_real_ else min(Re(eig_vals), na.rm=TRUE)
  eig_max <- if (all(is.na(eig_vals))) NA_real_ else max(Re(eig_vals), na.rm=TRUE)

  diagnostics[[as.character(window_size)]] <- list(
    window_size = window_size,
    n_snps = n_snps_clean,
    n_groups = n_groups,
    pool_members_sizes = purrr::map_int(pool_members, length),
    cov_dim = dim(cov_pool),
    kappa = kappa_val,
    eig_min = eig_min,
    eig_max = eig_max,
    diag_sum = diag_sum,
    off_sum = off_sum,
    na_dropped = na_dropped,
    sample_range = sample_range,
    founder_cor_offdiag_max = founder_cor_offdiag_max,
    founder_cor_offdiag_min = founder_cor_offdiag_min,
    A_pool_rank = rank_est,
    A_pool_cor_offdiag_max = ap_cor_offdiag_max,
    XtX_eig_min = XtX_eig_min,
    XtX_eig_max = XtX_eig_max,
    XtX_kappa = XtX_kappa,
    sigma2 = sigma2
  )
}

cat("Diagnostics by window size (bp):\n")
purrr::iwalk(diagnostics, function(d, nm){
  cat(sprintf("%6s: n_snps=%-5d n_groups=%-2d cov=%sx%s kappa=%.3e eig=[%.3e, %.3e] diag=%.3e off=%.3e | NA_drop=%d samp_rng=[%.3f, %.3f] f_cor_max=%.3f A_cor_max=%.3f A_rank=%s XtX_kappa=%.3e XtX_eig=[%.3e, %.3e] sigma2=%.3e\n",
              nm, d$n_snps, d$n_groups, d$cov_dim[1], d$cov_dim[2], d$kappa, d$eig_min, d$eig_max, d$diag_sum, d$off_sum,
              d$na_dropped, d$sample_range[1], d$sample_range[2], d$founder_cor_offdiag_max, d$A_pool_cor_offdiag_max,
              ifelse(is.na(d$A_pool_rank), "NA", as.character(d$A_pool_rank)), d$XtX_kappa, d$XtX_eig_min, d$XtX_eig_max, d$sigma2))
})

cat("\nDone.\n")


