# Just the est_haps_var function - no workflow code

suppressPackageStartupMessages({
  library(tidyverse)
  library(limSolve)
  library(MASS)
})

# Advanced haplotype estimator with genomic distance-based windowing and progressive V matrix
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
                     E = E, F = F, G = diag(n_founders), H = matrix(0, n_founders, 1),
                     fulloutput = TRUE),
      error = function(e) NULL
    )

    if (is.null(res)) {
      if (verbose >= 2) cat("    LSEI failed; continue\n")
      next
    }

    # Store design for true-cov computation
    last_A <- founder_matrix_clean
    last_y <- sample_freqs_clean

    # Compute pooled covariance
    pooled_result <- pooled_cov(founder_matrix_clean, sample_freqs_clean, groups)
    if (is.null(pooled_result$cov) || any(is.na(pooled_result$cov))) {
      if (verbose >= 2) cat("    Pooled cov failed; continue\n")
      next
    }

    # Update V matrix with pooled covariances
    for (i in seq_along(pooled_result$members)) {
      members <- pooled_result$members[[i]]
      for (ii in members) {
        for (jj in members) {
          V[ii, jj] <- pooled_result$cov[i, i]
        }
      }
    }

    # Store result
    final_result <- list(
      Groups = groups,
      Haps = setNames(as.numeric(res$X), founders),
      Err = V,
      Names = founders
    )

    if (verbose >= 2) {
      cat(sprintf("    Groups: %s\n", paste(group_strings, collapse = ", ")))
      cat(sprintf("    Haps: %s\n", paste(sprintf("%.3f", res$X), collapse = ", ")))
      print_V_compact(V)
    }
  }

  if (is.null(final_result)) {
    # Fallback: return minimal result
    final_result <- list(
      Groups = rep(1, length(founders)),
      Haps = setNames(rep(1/length(founders), length(founders)), founders),
      Err = diag(length(founders)),
      Names = founders
    )
  }

  return(final_result)
}
