#!/usr/bin/env Rscript

# =============================================================================
# WIDE FORMAT HAPLOTYPE ESTIMATOR WITH PROGRESSIVE ERROR MATRIX
# =============================================================================
# 
# This function takes wide format data and real genomic positions
# Uses the same progressive error matrix algorithm we developed
# Works with position-centered windows like production code

est_haps_wide <- function(pos, sample_name, df3_wide, founders, h_cutoff,
                         method = "adaptive",
                         window_size_bp = NULL,
                         chr = "chr2R",
                         verbose = 0) {
  
  # This is the progressive error matrix version adapted for wide format data
  # and real genomic positions (not simulated prefix-based windows)
  
  if (verbose >= 2) {
    cat(sprintf("=== WIDE FORMAT MODE: h_cutoff=%g, pos=%s, sample=%s ===\n", 
                h_cutoff, format(pos, big.mark=","), sample_name))
  }
  
  # Use production window sizes (in bp) and position-centered windows
  window_sizes <- c(10000, 20000, 50000, 100000, 200000, 500000)
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
  
  # Ensure POS-ordered wide input
  df3_wide <- df3_wide %>% dplyr::arrange(POS)
  
  for (window_size in window_sizes) {
    # Position-centered window (like production)
    window_start <- max(1, pos - window_size/2)
    window_end <- pos + window_size/2
    
    if (verbose >= 2) {
      cat(sprintf("  Window %s: %s", 
                  ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp")),
                  ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp"))))
    }
    
    # Subset wide data for this window
    window_data <- df3_wide %>%
      dplyr::filter(POS >= window_start & POS <= window_end)
    
    if (nrow(window_data) == 0) next
    
    # Check if we have the required columns
    required_cols <- c(founders, sample_name)
    if (!all(required_cols %in% names(window_data))) next
    
    # Extract founder matrix and sample frequencies directly from wide format
    founder_matrix <- window_data %>% 
      dplyr::select(dplyr::all_of(founders)) %>% 
      as.matrix()
    sample_freqs <- window_data %>% 
      dplyr::pull(!!sample_name)
    
    complete_rows <- stats::complete.cases(founder_matrix) & !is.na(sample_freqs)
    founder_matrix_clean <- founder_matrix[complete_rows, , drop = FALSE]
    sample_freqs_clean <- sample_freqs[complete_rows]
    
    if (nrow(founder_matrix_clean) < 10) next
    
    if (verbose >= 2) {
      cat(sprintf(" - %d SNPs, %d groups", nrow(founder_matrix_clean), n_groups))
    }
    
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
    
    # Check if clustering improved
    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat(" - ✗ No improvement\n")
      next  # No improvement, try larger window
    }
    
    previous_n_groups <- n_groups
    
    # Constraints accumulate like production
    n_founders <- ncol(founder_matrix_clean)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("    ✓ Added %d accumulated constraints from previous windows\n", nrow(accumulated_constraints)))
    }
    
    res <- tryCatch(
      limSolve::lsei(A = founder_matrix_clean, B = sample_freqs_clean,
                     E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat(" - ✗ LSEI failed\n")
      next
    }
    
    if (verbose >= 2) {
      cat(sprintf(" - ✓ LSEI success, %d groups", n_groups))
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
        cat(sprintf(" - %d constraints", built_ct))
      }
    } else {
      accumulated_constraints <- NULL
      accumulated_constraint_values <- NULL
      if (verbose >= 2) {
        built_ct <- 0
      }
    }
    
    final_result <- res
    final_n_groups <- n_groups
    groups_path <- c(groups_path, n_groups)
    
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
      # Approximate weights from current (last) res if same dimensionality
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
      if (verbose >= 2) cat(" - ✅ SUCCESS: All founders distinguished\n")
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
