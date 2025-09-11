#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
  library(MASS)
})

# simulate_founders
# Parent–child haplotype simulator generating an n_snps x n_founders (8) 0/1 matrix.
# Concept:
# - Draw K base (unrelated) founders using independent Bernoulli(0.5) SNPs.
# - The remaining founders are children derived from these bases.
# - Each child is assigned a target window size W ∈ {300, 750, 1500, 3000}.
# - For that child, apply per-150-SNP-block flips so that within window W,
#   the child differs from its parent by >20 SNPs on average.
#   flips_per_block = round(20 * 150 / W).
# - This creates predictable distinguishability timings: smaller W → earlier split.
#
# Parameters:
# - n_snps: total number of SNPs (e.g., 3000)
# - n_founders: number of founders (e.g., 8)
# Returns:
# - Integer matrix of shape [n_snps, n_founders] with entries in {0,1}.
simulate_founders <- function(n_snps, n_founders) {
  block <- 150L
  n_blocks <- n_snps %/% block
  # Number of unrelated bases (Poisson around 5.5), clamp to [3, 7] and ≤ n_founders
  K <- max(3L, min(7L, rpois(1, lambda = 5.5)))
  K <- min(K, n_founders)
  # Create K base templates
  bases <- replicate(K, rbinom(n_snps, 1, 0.5), simplify = FALSE)
  # Determine number of children and assign each to a parent among bases
  num_children <- n_founders - K
  parents <- if (num_children > 0) sample(seq_len(K), size = num_children, replace = TRUE) else integer(0)
  # Initialize founder matrix and place bases
  A <- matrix(0L, nrow = n_snps, ncol = n_founders)
  for (i in seq_len(K)) A[, i] <- bases[[i]]
  # Child mutation schedule by target window
  target_windows <- c(300L, 750L, 1500L, 3000L)
  calc_flips_per_block <- function(target_w) as.integer(round(20 * 150 / target_w))
  mutate_child <- function(parent_vec, flips_per_block) {
    out <- parent_vec
    if (flips_per_block <= 0) return(out)
    for (b in seq_len(n_blocks)) {
      s <- (b - 1L) * block + 1L; e <- min(b * block, n_snps)
      idx <- s:e
      k <- min(flips_per_block, length(idx))
      if (k > 0) {
        flip <- sample(idx, size = k)
        out[flip] <- 1L - out[flip]
      }
    }
    out
  }
  if (num_children > 0) {
    for (j in seq_len(num_children)) {
      parent_idx <- parents[j]
      target_w <- sample(target_windows, 1)
      flips_per_block <- calc_flips_per_block(target_w)
      A[, K + j] <- mutate_child(A[, parent_idx], flips_per_block)
    }
  }
  # Safety: if fewer columns filled, add random bases
  if (K + num_children < n_founders) {
    for (c in (K + num_children + 1):n_founders) A[, c] <- rbinom(n_snps, 1, 0.5)
  }
  A
}

# Do NOT source production here to avoid triggering its non-interactive main.
# Provide a stub alias; this workbench uses simulated prefix mode (pos = -99)
# and does not require calling production.
estimate_haplotypes_list_format_prod <- function(...) {
  stop("Production estimator not available in workbench; use pos = -99 to run simulated prefix mode.")
}

# Batch collector returning a tibble with metrics (no printing)
run_batch_df <- function(n_runs = 12, n_snps = 3000L, h_cutoff = 4) {
  out <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    set.seed(1000 + r)
    n_founders <- 8
    founders <- paste0("F", seq_len(n_founders))
    sample_name <- "S1"
    POS <- seq_len(n_snps)

    A <- matrix(0, nrow = n_snps, ncol = n_founders)
    block <- 150L; n_blocks <- n_snps %/% block
    flip_blocks <- function(col, k) {
      for (b in seq_len(n_blocks)) {
        s <- (b - 1L) * block + 1L; e <- min(b * block, n_snps)
        idx <- s:e
        if (k > 0) { flip <- sample(idx, size = min(k, length(idx))); col[flip] <- 1 - col[flip] }
      }
      col
    }
    set.seed(1100 + r)
    B1 <- rbinom(n_snps, 1, 0.5); B2 <- rbinom(n_snps, 1, 0.5); B3 <- rbinom(n_snps, 1, 0.5)
    bases <- list(B1 = B1, B2 = B2, B3 = B3)
    repeat {
      base_assign <- sample(c("B1", "B2", "B3"), size = n_founders, replace = TRUE, prob = c(0.34, 0.33, 0.33))
      tab <- table(base_assign)
      if (all(tab <= 2)) break
    }
    if (length(which(table(base_assign) == 2)) == 0) base_assign[2] <- base_assign[1]
    flip_opts <- c(4, 6, 8, 10, 12)
    flips <- integer(n_founders)
    for (b in names(bases)) {
      idx <- which(base_assign == b)
      if (length(idx) == 1) {
        flips[idx] <- sample(flip_opts, 1)
      } else if (length(idx) == 2) {
        fpair <- sample(flip_opts, 2, replace = FALSE)
        while (abs(diff(fpair)) < 4) fpair <- sample(flip_opts, 2, replace = FALSE)
        flips[idx] <- fpair
      }
    }
    for (f in seq_len(n_founders)) A[, f] <- flip_blocks(bases[[ base_assign[f] ]], flips[f])
    colnames(A) <- founders

    set.seed(2000 + r)
    w <- runif(n_founders); w <- w / sum(w)
    noise_sd <- ifelse(r %% 3 == 0, 0.03, 0.02)
    samp <- as.numeric(A %*% w) + rnorm(n_snps, 0, noise_sd)
    samp[samp < 0] <- 0; samp[samp > 1] <- 1

    df_founders <- tibble::as_tibble(A) %>% dplyr::mutate(POS = POS) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(founders), names_to = "name", values_to = "freq")
    df_sample <- tibble::tibble(POS = POS, name = sample_name, freq = samp)
    df3 <- dplyr::bind_rows(df_founders, df_sample)

    res <- estimate_haplotypes_list_format_sim(
      pos = -99, sample_name = sample_name, df3 = df3, founders = founders,
      h_cutoff = h_cutoff, method = "adaptive", verbose = 0
    )

    hap_err <- max(abs(as.numeric(res$Haps) - w)) * 100
    ek <- tryCatch(kappa(res$Err), error=function(e) NA_real_)
    n_groups_final <- length(unique(res$Groups))
    out[[r]] <- tibble::tibble(run = r, hap_err = hap_err, kappa = ek, Ng = n_groups_final)
  }
  dplyr::bind_rows(out)
}
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
  groups_path <- integer(0)
  groups_fmt_path <- character(0)
  groups_win_path <- integer(0)

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  # Progressive V matrix (8x8), pending covariances, and resolved tracker
  V <- matrix(NA_real_, length(founders), length(founders))
  rownames(V) <- founders; colnames(V) <- founders
  pending <- list()
  resolved <- setNames(rep(FALSE, length(founders)), founders)

  # Track last (largest successful) unpooled design for true-cov computation
  last_A <- NULL
  last_y <- NULL

  # Helpers to compute pooled covariance at a window and to print V compactly
  pooled_cov <- function(A_full, y_full, groups_vec){
    gid <- sort(unique(groups_vec))
    k <- length(gid)
    A_pool <- matrix(0, nrow(A_full), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups_vec == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
    }
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar)) {
      return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
    }
    XtX <- crossprod(A_pool)
    xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) ginv(XtX) %*% crossprod(A_pool, y_full))
    r <- y_full - as.numeric(A_pool %*% xhat)
    sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
    list(cov = sigma2 * tryCatch(solve(XtX), error=function(e) ginv(XtX)), members=members, w=as.numeric(xhat))
  }

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
    m2 <- as.integer(round(mant * 10))  # two-digit mantissa (approx)
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

  # Ensure POS-ordered tidy input (POS, name, freq)
  df3 <- df3 %>% dplyr::arrange(POS)

  for (window_size in window_sizes) {
    window_data <- df3 %>%
      dplyr::filter(POS >= 1, POS <= window_size, name %in% c(founders, sample_name))
    if (nrow(window_data) == 0) next

    wide_data <- window_data %>%
      dplyr::select(POS, name, freq) %>%
      tidyr::pivot_wider(names_from = name, values_from = freq)

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
    if (verbose >= 2) {
      comp <- split(founders, groups)
      ordered_gids <- sort(as.integer(names(comp)))
      group_strings <- vapply(ordered_gids, function(gid){
        paste(comp[[as.character(gid)]], collapse="+")
      }, character(1))
      carried_ct <- if (!is.null(accumulated_constraints)) nrow(accumulated_constraints) else 0
    }

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
      if (verbose >= 2) cat(sprintf("    Carried over %d constraints\n", nrow(accumulated_constraints)))
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
      if (verbose >= 2) {
        built_ct <- nrow(current_constraints)
        # Map values to ordered group ids
        ordered_gids <- sort(unique(groups))
        group_values <- vapply(ordered_gids, function(gid){
          sum(res$X[which(groups == gid)])
        }, numeric(1))
      }
    } else {
      accumulated_constraints <- NULL
      accumulated_constraint_values <- NULL
      if (verbose >= 2) {
        built_ct <- 0
        group_values <- numeric(0)
      }
    }

    final_result <- res
    final_n_groups <- n_groups
    groups_path <- c(groups_path, n_groups)
    # compact groups format like [123]4[56]78 based on founders order
    fmt <- {
      grp_ids <- sort(unique(groups))
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
    groups_fmt_path <- c(groups_fmt_path, fmt)
    groups_win_path <- c(groups_win_path, window_size)
    if (verbose >= 2) {
      # Print compact, aligned 3-line block (monospace-friendly)
      cat(sprintf("window_snp=%-5d  n_snps=%-5d  n_groups=%-2d  carried=%-2d  built=%-2d\n",
                  window_size, nrow(founder_matrix_clean), n_groups, carried_ct, built_ct))
      if (length(group_strings)) {
        # Column widths match label lengths exactly for alignment
        col_w <- nchar(group_strings)
        # Render group labels with padding
        grp_fmt <- vapply(seq_along(group_strings), function(i){
          sprintf("%-*s", col_w[i], group_strings[i])
        }, character(1))
        cat(paste(grp_fmt, collapse=" | "), "\n", sep="")
      }
      if (length(group_values)) {
        # Convert to percent, 0 decimals
        vals_pct <- round(group_values * 100, 0)
        # Match widths of labels exactly
        if (!exists("col_w")) col_w <- rep(1, length(vals_pct))
        val_fmt <- vapply(seq_along(vals_pct), function(i){
          sprintf("%*d", col_w[i], as.integer(vals_pct[i]))
        }, character(1))
        cat(paste(val_fmt, collapse=" | "), "\n", sep="")
      }
    }

    # Update progressive V using pooled model for this window
    pc <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
    cov_pool <- pc$cov; pool_members <- pc$members
    # mark newly resolved founders
    for (ii in seq_along(pool_members)){
      mem <- pool_members[[ii]]
      if (length(mem)==1L) resolved[founders[mem]] <- TRUE
    }
    # write cov for resolved-resolved; store pending for det-vs-pool and pool-vs-pool
    for (i_idx in seq_along(pool_members)){
      mi <- pool_members[[i_idx]]
      for (j_idx in seq_along(pool_members)){
        mj <- pool_members[[j_idx]]
        if (length(mi)==1L && length(mj)==1L){
          fi <- founders[mi]; fj <- founders[mj]
          if (is.na(V[fi,fj])) V[fi,fj] <- cov_pool[i_idx, j_idx]
          if (is.na(V[fj,fi])) V[fj,fi] <- cov_pool[j_idx, i_idx]
        } else if (length(mi)==1L && length(mj)>1L){
          pending <- append(pending, list(list(type="det_vs_pool", det=founders[mi], pool=j_idx, cov=cov_pool[i_idx, j_idx], mem=founders[mj])))
        } else if (length(mi)>1L && length(mj)>1L && i_idx<=j_idx){
          pending <- append(pending, list(list(type="pool_vs_pool", pool_a=i_idx, pool_b=j_idx, cov=cov_pool[i_idx, j_idx], mem_a=founders[mi], mem_b=founders[mj])))
        }
      }
    }

    # If all resolved now, distribute pending covariances using current weights
    if (all(resolved)){
      w_now <- setNames(rep(NA_real_, length(founders)), founders)
      # Approximate weights from current (last) res if same dimensionality, else use group_values map
      if (length(final_result$X) == length(founders)){
        w_now <- setNames(as.numeric(final_result$X), founders)
      }
      for (rec in pending){
        if (rec$type == "det_vs_pool"){
          den <- sum(w_now[rec$mem]); if (!is.finite(den) || den<=0) next
          for (f in rec$mem){
            val <- rec$cov * (w_now[f]/den)
            if (is.na(V[rec$det, f])) V[rec$det, f] <- val
            if (is.na(V[f, rec$det])) V[f, rec$det] <- val
          }
        } else if (rec$type == "pool_vs_pool"){
          den_a <- sum(w_now[rec$mem_a]); den_b <- sum(w_now[rec$mem_b]); if (den_a<=0 || den_b<=0) next
          for (fa in rec$mem_a){
            for (fb in rec$mem_b){
              val <- rec$cov * (w_now[fa]/den_a) * (w_now[fb]/den_b)
              if (is.na(V[fa, fb])) V[fa, fb] <- val
              if (is.na(V[fb, fa])) V[fb, fa] <- val
            }
          }
        }
      }
      pending <- list()
    }

    if (verbose >= 2){
      print_V_compact(V)
    }
    if (n_groups == length(founders)) {
      # Save unpooled design at the final successful window
      last_A <- founder_matrix_clean
      last_y <- sample_freqs_clean
      if (verbose >= 2) cat("    All founders distinguished; stop\n")
      break
    }
  }

  # Return EXACT same structure as production
  if (!is.null(final_result)) {
    founder_frequencies <- final_result$X; names(founder_frequencies) <- founders
    # Compute "true" covariance using largest window residual sigma^2 and unpooled design
    Cov_true <- NULL
    if (!is.null(last_A) && !is.null(last_y)) {
      XtX <- crossprod(last_A)
      Xinv <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
      r <- last_y - as.numeric(last_A %*% founder_frequencies)
      p <- ncol(last_A); n <- nrow(last_A)
      sigma2_hat <- sum(r^2) / max(1, n - p)
      Cov_true <- sigma2_hat * Xinv
    }
    # Prefer progressive V if sufficiently filled; otherwise fallback
    if (sum(is.na(V)) < length(V)) {
      for (d in seq_len(nrow(V))) if (is.na(V[d, d])) V[d, d] <- 1e-8
      error_matrix <- V
    } else if (!is.null(Cov_true)) {
      error_matrix <- Cov_true
    } else if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
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

  res_out <- list(Groups=groups, Haps=founder_frequencies, Err=error_matrix, Names=founders)
  # Attach diagnostics for wrapper comparison
  attr(res_out, "V_progressive") <- V
  if (exists("Cov_true")) attr(res_out, "Cov_true") <- Cov_true else attr(res_out, "Cov_true") <- NULL
  attr(res_out, "groups_path") <- groups_path
  attr(res_out, "groups_fmt_path") <- groups_fmt_path
  attr(res_out, "groups_win_path") <- groups_win_path
  res_out
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
    # Local compact printer for comparison output
    fmt_cell_signed_local <- function(x, diag_cell){
      if (is.na(x)) return("  NA ")
      if (x == 0) { signc <- if (diag_cell) " " else "+"; return(paste0(signc, "  0 ")) }
      signc <- if (diag_cell) " " else if (x < 0) "-" else "+"
      ax <- abs(x); expo <- floor(log10(ax)); mant <- ax / (10^expo)
      m2 <- as.integer(round(mant * 10)); if (expo < -9) expo <- -9; if (expo > 9) expo <- 9
      paste0(signc, sprintf("%02d%+1d", m2, as.integer(expo)))
    }
    print_matrix_compact <- function(M){
      for (i in seq_len(nrow(M))){
        row <- character(ncol(M))
        for (j in seq_len(ncol(M))) row[j] <- fmt_cell_signed_local(M[i,j], diag_cell=(i==j))
        cat(paste(row, collapse=" "), "\n", sep="")
      }
    }
    cat("Groups:", paste(res$Groups, collapse=","), "\n")
    # Pretty print true vs estimated hap frequencies as aligned table (percent, 0 decimals)
    cat("True vs Estimated Haps (%)\n")
    labs <- founders
    vals_true <- round(w * 100)
    vals_est  <- round(as.numeric(res$Haps) * 100)
    col_w <- pmax(nchar(labs), 2)
    lab_line <- paste(vapply(seq_along(labs), function(i) sprintf("%-*s", col_w[i], labs[i]), character(1)), collapse=" | ")
    tru_line <- paste(vapply(seq_along(vals_true), function(i) sprintf("%*d", col_w[i], vals_true[i]), character(1)), collapse=" | ")
    est_line <- paste(vapply(seq_along(vals_est),  function(i) sprintf("%*d", col_w[i], vals_est[i]),  character(1)), collapse=" | ")
    cat(lab_line, "\n", sep="")
    cat(tru_line, "\n", sep="")
    cat(est_line, "\n", sep="")
    if (is.matrix(res$Err)) cat("Err kappa:", kappa(res$Err), "\n")
    V_est <- attr(res, "V_progressive")
    Cov_true <- attr(res, "Cov_true")
    if (!is.null(V_est)) { cat("Estimated V (signed, sci):\n"); print_matrix_compact(V_est) }
    if (!is.null(Cov_true)) { cat("Cov_true (signed, sci):\n"); print_matrix_compact(Cov_true) }
    if (!is.null(V_est) && !is.null(Cov_true)) {
      diffs <- as.numeric(V_est - Cov_true)
      diffs <- diffs[is.finite(diffs)]
      cat(sprintf("Diff summary | max abs: %.3e  median abs: %.3e\n", max(abs(diffs)), median(abs(diffs))))
      cat(sprintf("Cond(V_est)=%.2e  Cond(Cov_true)=%.2e\n", kappa(V_est), kappa(Cov_true)))
    }
  }

  invisible(res)
}

# Batch runner: minimal tabular output over multiple random runs
run_batch <- function(n_runs = 12, n_snps = 3000L, h_cutoff = 4) {
  # header
  cat(sprintf("%4s  %6s  %8s  %3s  %s\n", "run", "hap%", "kappa", "Ng", "groups progression (win:groups)"))
  for (r in seq_len(n_runs)) {
    # vary seeds for diversity
    set.seed(1000 + r)
    n_founders <- 8
    founders <- paste0("F", seq_len(n_founders))
    sample_name <- "S1"
    POS <- seq_len(n_snps)

    # Parent–child simulator per spec (documented above)
    set.seed(1100 + r)
    A <- simulate_founders(n_snps, n_founders)
    colnames(A) <- founders
    colnames(A) <- founders

    set.seed(2000 + r)
    w <- runif(n_founders); w <- w / sum(w)
    noise_sd <- ifelse(r %% 3 == 0, 0.03, 0.02)
    samp <- as.numeric(A %*% w) + rnorm(n_snps, 0, noise_sd)
    samp[samp < 0] <- 0; samp[samp > 1] <- 1

    df_founders <- tibble::as_tibble(A) %>% dplyr::mutate(POS = POS) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(founders), names_to = "name", values_to = "freq")
    df_sample <- tibble::tibble(POS = POS, name = sample_name, freq = samp)
    df3 <- dplyr::bind_rows(df_founders, df_sample)

    res <- estimate_haplotypes_list_format_sim(
      pos = -99, sample_name = sample_name, df3 = df3, founders = founders,
      h_cutoff = h_cutoff, method = "adaptive", verbose = 0
    )

    # summarize groups as compact strings
    gpf <- attr(res, "groups_fmt_path"); gw <- attr(res, "groups_win_path")
    if (!is.null(gw) && length(gw) == length(gpf)) {
      gp_str <- paste(sprintf("%d:%s", gw, gpf), collapse = " -> ")
    } else {
      gp_str <- paste(gpf, collapse = " -> ")
    }
    hap_err <- max(abs(as.numeric(res$Haps) - w)) * 100
    ek <- tryCatch(kappa(res$Err), error=function(e) NA_real_)
    # Format kappa: plain one-decimal if <1000; scientific if >=1000; Inf/NA as "Inf"
    kappa_str <- if (is.na(ek) || !is.finite(ek)) {
      "Inf"
    } else if (ek < 1000) {
      sprintf("%.1f", ek)
    } else {
      sprintf("%.1e", ek)
    }
    n_groups_final <- length(unique(res$Groups))
    cat(sprintf("%4d  %6.1f  %8s  %3d  %s\n", r, hap_err, kappa_str, n_groups_final, gp_str))
  }
}

# Replay a single run index with verbose diagnostics to investigate issues
run_one_verbose <- function(run_index = 3, n_snps = 3000L, h_cutoff = 4, verbose = 2) {
  set.seed(1000 + run_index)
  n_founders <- 8
  founders <- paste0("F", seq_len(n_founders))
  sample_name <- "S1"
  POS <- seq_len(n_snps)

  A <- matrix(0, nrow = n_snps, ncol = n_founders)
  block <- 150L; n_blocks <- n_snps %/% block
  flip_blocks <- function(col, k) {
    for (b in seq_len(n_blocks)) {
      s <- (b - 1L) * block + 1L; e <- min(b * block, n_snps)
      idx <- s:e
      if (k > 0) { flip <- sample(idx, size = min(k, length(idx))); col[flip] <- 1 - col[flip] }
    }
    col
  }
  # Parent–child simulator for replay (documented above)
  set.seed(1100 + run_index)
  A <- simulate_founders(n_snps, n_founders)
  colnames(A) <- founders

  set.seed(2000 + run_index)
  w <- runif(n_founders); w <- w / sum(w)
  noise_sd <- ifelse(run_index %% 3 == 0, 0.03, 0.02)
  samp <- as.numeric(A %*% w) + rnorm(n_snps, 0, noise_sd)
  samp[samp < 0] <- 0; samp[samp > 1] <- 1

  df_founders <- tibble::as_tibble(A) %>% dplyr::mutate(POS = POS) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(founders), names_to = "name", values_to = "freq")
  df_sample <- tibble::tibble(POS = POS, name = sample_name, freq = samp)
  df3 <- dplyr::bind_rows(df_founders, df_sample)

  res <- estimate_haplotypes_list_format_sim(
    pos = -99, sample_name = sample_name, df3 = df3, founders = founders,
    h_cutoff = h_cutoff, method = "adaptive", verbose = verbose
  )
  res
}


# Convenience runner: print 100-run table then a concise summary
run_batch_100_with_summary <- function(h_cutoff = 4) {
  run_batch(100, h_cutoff = h_cutoff)
  df <- run_batch_df(100, h_cutoff = h_cutoff)
  df <- df %>% dplyr::mutate(kappa_inf = !is.finite(kappa))
  n <- nrow(df)
  inf <- sum(df$kappa_inf)
  ok <- n - inf
  all_inf_bad <- all(ifelse(df$kappa_inf, df$Ng != 8, TRUE))
  min_ng_inf <- if (inf > 0) min(df$Ng[df$kappa_inf]) else NA_integer_
  max_ng_inf <- if (inf > 0) max(df$Ng[df$kappa_inf]) else NA_integer_
  cat("\nSummary (100 runs):\n", sep="")
  cat(sprintf("  total=%d  kappa_Inf=%d  kappa_ok=%d  all_Inf_have_Ng!=8=%s\n",
              n, inf, ok, ifelse(all_inf_bad, "TRUE", "FALSE")))
  if (inf > 0) cat(sprintf("  Ng among Inf: min=%d max=%d\n", min_ng_inf, max_ng_inf))
}

# Enhanced batch collector: returns tibble with groups progression as text column
run_batch_df_enhanced <- function(n_runs = 12, n_snps = 3000L, h_cutoff = 4) {
  out <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    set.seed(1000 + r)
    n_founders <- 8
    founders <- paste0("F", seq_len(n_founders))
    sample_name <- "S1"
    POS <- seq_len(n_snps)

    set.seed(1100 + r)
    A <- simulate_founders(n_snps, n_founders)
    colnames(A) <- founders

    set.seed(2000 + r)
    w <- runif(n_founders); w <- w / sum(w)
    noise_sd <- ifelse(r %% 3 == 0, 0.03, 0.02)
    samp <- as.numeric(A %*% w) + rnorm(n_snps, 0, noise_sd)
    samp[samp < 0] <- 0; samp[samp > 1] <- 1

    df_founders <- tibble::as_tibble(A) %>% dplyr::mutate(POS = POS) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(founders), names_to = "name", values_to = "freq")
    df_sample <- tibble::tibble(POS = POS, name = sample_name, freq = samp)
    df3 <- dplyr::bind_rows(df_founders, df_sample)

    res <- estimate_haplotypes_list_format_sim(
      pos = -99, sample_name = sample_name, df3 = df3, founders = founders,
      h_cutoff = h_cutoff, method = "adaptive", verbose = 0
    )

    hap_err <- max(abs(as.numeric(res$Haps) - w)) * 100
    ek <- tryCatch(kappa(res$Err), error=function(e) NA_real_)
    n_groups_final <- length(unique(res$Groups))
    
    # Build groups progression string
    gpf <- attr(res, "groups_fmt_path"); gw <- attr(res, "groups_win_path")
    groups_progression <- if (!is.null(gw) && length(gw) == length(gpf)) {
      paste(sprintf("%d:%s", gw, gpf), collapse = " -> ")
    } else {
      paste(gpf, collapse = " -> ")
    }
    
    out[[r]] <- tibble::tibble(
      run = r, 
      hap_err = hap_err, 
      kappa = ek, 
      Ng = n_groups_final,
      groups_progression = groups_progression
    )
  }
  dplyr::bind_rows(out)
}

# Summary functions for batch results
summarize_batch_results <- function(df) {
  df <- df %>% dplyr::mutate(
    kappa_inf = !is.finite(kappa),
    converged = Ng == 8,
    kappa_log10 = ifelse(kappa_inf, NA_real_, log10(kappa))
  )
  
  # Basic counts
  n <- nrow(df)
  inf <- sum(df$kappa_inf)
  ok <- n - inf
  converged <- sum(df$converged)
  not_converged <- n - converged
  
  # Validation: Inf only when not converged
  all_inf_bad <- all(ifelse(df$kappa_inf, !df$converged, TRUE))
  
  # Hap error stats
  hap_stats <- df %>% dplyr::summarise(
    hap_mean = mean(hap_err, na.rm = TRUE),
    hap_median = median(hap_err, na.rm = TRUE),
    hap_sd = sd(hap_err, na.rm = TRUE),
    hap_min = min(hap_err, na.rm = TRUE),
    hap_max = max(hap_err, na.rm = TRUE)
  )
  
  # Kappa stats (finite only)
  kappa_stats <- df %>% 
    dplyr::filter(!kappa_inf) %>% 
    dplyr::summarise(
      kappa_mean = mean(kappa_log10, na.rm = TRUE),
      kappa_median = median(kappa_log10, na.rm = TRUE),
      kappa_sd = sd(kappa_log10, na.rm = TRUE),
      kappa_min = min(kappa_log10, na.rm = TRUE),
      kappa_max = max(kappa_log10, na.rm = TRUE)
    )
  
  # Groups progression analysis
  progression_stats <- df %>% 
    dplyr::mutate(
      n_transitions = stringr::str_count(groups_progression, " -> "),
      reaches_8 = stringr::str_detect(groups_progression, "12345678$")
    ) %>% 
    dplyr::summarise(
      mean_transitions = mean(n_transitions, na.rm = TRUE),
      pct_reach_8 = mean(reaches_8, na.rm = TRUE) * 100
    )
  
  list(
    counts = list(
      total = n, 
      kappa_inf = inf, 
      kappa_finite = ok,
      converged = converged,
      not_converged = not_converged,
      all_inf_bad = all_inf_bad
    ),
    hap_error = hap_stats,
    kappa_log10 = kappa_stats,
    progression = progression_stats
  )
}

print_batch_summary <- function(df) {
  s <- summarize_batch_results(df)
  
  cat("=== Batch Results Summary ===\n")
  cat(sprintf("Runs: %d total\n", s$counts$total))
  cat(sprintf("Convergence: %d converged (Ng=8), %d not converged\n", 
              s$counts$converged, s$counts$not_converged))
  cat(sprintf("Kappa: %d finite, %d infinite\n", 
              s$counts$kappa_finite, s$counts$kappa_inf))
  cat(sprintf("Inf kappa only when not converged: %s\n", 
              ifelse(s$counts$all_inf_bad, "TRUE", "FALSE")))
  
  cat("\nHap Error (%):\n")
  cat(sprintf("  mean=%.2f  median=%.2f  sd=%.2f  range=[%.2f,%.2f]\n",
              s$hap_error$hap_mean, s$hap_error$hap_median, s$hap_error$hap_sd,
              s$hap_error$hap_min, s$hap_error$hap_max))
  
  if (s$counts$kappa_finite > 0) {
    cat("\nLog10(Kappa) (finite only):\n")
    cat(sprintf("  mean=%.2f  median=%.2f  sd=%.2f  range=[%.2f,%.2f]\n",
                s$kappa_log10$kappa_mean, s$kappa_log10$kappa_median, s$kappa_log10$kappa_sd,
                s$kappa_log10$kappa_min, s$kappa_log10$kappa_max))
  }
  
  cat("\nGroups Progression:\n")
  cat(sprintf("  mean transitions=%.1f  %% reach 8 groups=%.1f\n",
              s$progression$mean_transitions, s$progression$pct_reach_8))
}

# Wrapper: run 100 simulations, show data frame, then summary
run_100_with_dataframe <- function(h_cutoff = 4) {
  cat("Running 100 simulations...\n")
  df <- run_batch_df_enhanced(100, h_cutoff = h_cutoff)
  
  cat("\n=== Data Frame (first 20 rows) ===\n")
  print(df[1:20, ])
  if (nrow(df) > 20) {
    cat(sprintf("... and %d more rows\n", nrow(df) - 20))
  }
  
  cat("\n")
  print_batch_summary(df)
  
  invisible(df)
}

