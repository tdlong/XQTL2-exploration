#!/usr/bin/env Rscript

# Debug wrapper - loads hunk data and calls exact production est_haps_var function
# with maximum verbosity to see what's happening

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(limSolve)
})

# Copy the EXACT est_haps_var function from production
est_haps_var <- function(testing_position, sample_name, df3, founders, h_cutoff,
                               method = "adaptive",
                               window_size_bp = NULL,
                               chr = "chr2R",
                               verbose = 0) {
  
  if (verbose >= 2) {
    cat(sprintf("=== ADVANCED HAPLOTYPE ESTIMATION: h_cutoff=%g, pos=%s, sample=%s ===\n", 
                h_cutoff, format(testing_position, big.mark=","), sample_name))
  }

  # Use EHLF window sizes (genomic distance-based)
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
    
    # Handle edge case where k=0 or k=1
    if (k == 0) {
      return(list(cov=matrix(NA, 0, 0), members=list(), w=numeric(0)))
    }
    if (k == 1) {
      return(list(cov=matrix(1, 1, 1), members=list(1), w=1))
    }
    
    A_pool <- matrix(0, nrow(A_full), k)
    members <- vector("list", k)
    for (ii in seq_along(gid)){
      mem <- which(groups_vec == gid[ii])
      members[[ii]] <- mem
      A_pool[, ii] <- if (length(mem) == 1L) A_full[, mem] else rowMeans(A_full[, mem, drop=FALSE])
    }
    
    E <- matrix(1, 1, k); F <- 1
    fit <- tryCatch(lsei(A=A_pool, B=y_full, E=E, F=F, G=diag(k), H=matrix(0, k, 1), fulloutput=TRUE), error=function(e) NULL)
    if (!is.null(fit) && !is.null(fit$covar) && is.matrix(fit$covar) && nrow(fit$covar) == k && ncol(fit$covar) == k) {
      return(list(cov=fit$covar, members=members, w=as.numeric(fit$X)))
    }
    
    # Fallback: compute covariance manually
    XtX <- crossprod(A_pool)
    xhat <- tryCatch(solve(XtX, crossprod(A_pool, y_full)), error=function(e) MASS::ginv(XtX) %*% crossprod(A_pool, y_full))
    r <- y_full - as.numeric(A_pool %*% xhat)
    sigma2 <- sum(r^2) / max(1, nrow(A_pool) - ncol(A_pool))
    
    cov_matrix <- tryCatch(solve(XtX), error=function(e) MASS::ginv(XtX))
    
    # Ensure we return a proper k×k matrix
    if (!is.matrix(cov_matrix) || nrow(cov_matrix) != k || ncol(cov_matrix) != k) {
      # Fallback: return identity matrix
      cov_matrix <- diag(k)
    }
    
    list(cov = sigma2 * cov_matrix, members=members, w=as.numeric(xhat))
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

  # Ensure POS-ordered wide input (POS, founder1, founder2, ..., sample1, sample2, ...)
  # df3 is now in WIDE format: POS, founder1, founder2, ..., foundern, sample1, sample2, ..., sampleM
 
  for (window_size in window_sizes) {
    # Calculate window boundaries based on genomic distance (like EHLF)
    window_start <- max(1, testing_position - window_size/2)
    window_end <- testing_position + window_size/2
    
    if (verbose >= 2) {
      cat(sprintf("  Window %s: %s", 
                  ifelse(window_size >= 1000, paste0(window_size/1000, "kb"), paste0(window_size, "bp")),
                  paste0("pos ", window_start, "-", window_end)))
    }
    
    # Filter to window positions (wide format - no pivoting needed!)
    window_data <- df3 %>%
      dplyr::filter(POS >= window_start & POS <= window_end)
    if (nrow(window_data) == 0) next

    # Check if all required samples are present
    if (!all(c(founders, sample_name) %in% names(window_data)) || nrow(window_data) < 10) next

    # Extract founder matrix and sample frequencies directly from wide format
    founder_matrix <- window_data %>% dplyr::select(dplyr::all_of(founders)) %>% as.matrix()
    sample_freqs <- window_data %>% dplyr::pull(!!sample_name)

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
    
    # Check for linear dependence before adding constraints
    if (!is.null(current_constraints)) {
      # Check if new constraints are linearly independent from existing ones
      if (!is.null(accumulated_constraints)) {
        # Combine existing and new constraints
        E_combined <- rbind(accumulated_constraints, current_constraints)
        
        # Check rank deficiency - if rank < nrow, there's linear dependence
        if (qr(E_combined)$rank < nrow(E_combined)) {
          if (verbose >= 2) {
            cat("    Warning: Linear dependence detected, skipping constraint addition\n")
            cat("    Existing constraints:", nrow(accumulated_constraints), "\n")
            cat("    New constraints:", nrow(current_constraints), "\n")
            cat("    Combined rank:", qr(E_combined)$rank, "vs rows:", nrow(E_combined), "\n")
          }
          # Skip adding these constraints to avoid numerical issues
          current_constraints <- NULL
          current_constraint_values <- NULL
        }
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
    
    # Skip V matrix update if only 1 group - no meaningful variance estimation possible
    if (n_groups == 1) {
      if (verbose >= 2) {
        cat("    Warning: Only 1 group, skipping V matrix update (no meaningful variance estimation)\n")
      }
      next  # Skip this window - no variance estimation possible with single group
    }
    
    # Debug: Check if cov_pool is a proper matrix
    if (!is.matrix(cov_pool)) {
      if (verbose >= 2) {
        cat("    Warning: cov_pool is not a matrix, skipping V matrix update\n")
        cat("    cov_pool class:", class(cov_pool), "\n")
        cat("    cov_pool length:", length(cov_pool), "\n")
      }
      next  # Skip this window if covariance is not computable
    }
    
    if (nrow(cov_pool) != length(pool_members) || ncol(cov_pool) != length(pool_members)) {
      if (verbose >= 2) {
        cat("    Warning: cov_pool dimensions don't match pool_members, skipping V matrix update\n")
        cat("    cov_pool dim:", dim(cov_pool), "pool_members length:", length(pool_members), "\n")
      }
      next  # Skip this window if dimensions don't match
    }
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

# Load hunk data
hunk_file <- "hunk_data_chr3R_19780000_Rep01_W_F_h4.rds"
if (!file.exists(hunk_file)) {
  stop("Hunk data file not found: ", hunk_file)
}

hunk_data <- readRDS(hunk_file)
cat("=== LOADING HUNK DATA ===\n")
cat("df3 dimensions:", nrow(hunk_data$df3), "x", ncol(hunk_data$df3), "\n")
cat("Args:", paste(names(hunk_data$args), collapse = ", "), "\n\n")

# Run with maximum verbosity
cat("=== RUNNING PRODUCTION ESTIMATOR WITH MAX VERBOSITY ===\n")
args_with_verbose <- hunk_data$args
args_with_verbose$verbose <- 3
result <- do.call(est_haps_var, c(list(df3 = hunk_data$df3), args_with_verbose))

cat("\n=== FINAL RESULT ===\n")
cat("Groups:", paste(result$Groups, collapse = ", "), "\n")
cat("Haps:", paste(sprintf("%.6f", result$Haps), collapse = ", "), "\n")
cat("Error matrix diagonal sum:", sum(diag(result$Err)), "\n")
cat("Error matrix diagonal:", paste(sprintf("%.2e", diag(result$Err)), collapse = ", "), "\n")

# Save result
saveRDS(result, "debug_result_chr3R_19780000_Rep01_W_F.rds")
cat("\nSaved debug result to: debug_result_chr3R_19780000_Rep01_W_F.rds\n")
