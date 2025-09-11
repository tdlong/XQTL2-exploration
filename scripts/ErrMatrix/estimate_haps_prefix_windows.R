#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limSolve)
})

##
## estimate_haplotypes_prefix
##
## Purpose:
## - Estimate haplotype proportions using the working adaptive LSEI approach,
##   but with prefix-based window growth (150, 300, 750, 1500, 3000 SNPs)
##   suitable for simulated data that is indexed linearly from 1..N.
## - Returns the SAME structure that production expects: list(Groups, Haps, Err, Names).
##
## How this maps back to production (position-centered windows):
## - Production centers windows around a genomic position and grows by bp lengths.
## - This prefix-based variant isolates ONLY the window-growth/solve logic.
## - To revert to the production calling style, use the adapter wrapper
##   `estimate_haplotypes_prefix_adapter()` below, which accepts the
##   production-like signature and extracts matrices accordingly.
## - The core returned structure and constraint accumulation behavior match production,
##   so swapping window selection (prefix → center) is trivial.
##
## Parameters:
## - founder_matrix: numeric matrix [num_snps × num_founders] of founder genotypes (0/1)
## - sample_freqs: numeric vector [num_snps] with sample allele frequencies in [0,1]
## - founder_names: character vector of length num_founders
## - h_cutoff: numeric, Euclidean cutoff passed to cutree
## - verbose: integer 0..2 (2 = detailed)
##
## Returns:
## - list(Groups, Haps, Err, Names)
##   - Groups: integer vector (cluster id per founder)
##   - Haps: named numeric vector of founder proportions
##   - Err: covariance matrix if provided by lsei (covar), otherwise NA matrix
##   - Names: character vector founder_names
##
## Example (simulated):
##   res <- estimate_haplotypes_prefix(A, y, founder_names, h_cutoff = 4, verbose = 2)
##

estimate_haplotypes_prefix <- function(
  founder_matrix,          # matrix [num_snps x num_founders] of 0/1 founder genotypes
  sample_freqs,            # numeric vector [num_snps] sample allele frequencies
  founder_names,           # character vector length num_founders
  h_cutoff = 4,            # clustering cutoff (Euclidean)
  verbose = 0              # 0-2
) {
  stopifnot(is.matrix(founder_matrix))
  stopifnot(length(sample_freqs) == nrow(founder_matrix))
  stopifnot(length(founder_names) == ncol(founder_matrix))

  window_sizes <- c(150, 300, 750, 1500, 3000)

  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founder_names))

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  for (window_size in window_sizes) {
    if (window_size > nrow(founder_matrix)) next

    A <- founder_matrix[seq_len(window_size), , drop = FALSE]
    y <- sample_freqs[seq_len(window_size)]

    # Clean NAs if any (shouldn't exist in simulated data)
    complete_rows <- stats::complete.cases(A) & !is.na(y)
    A <- A[complete_rows, , drop = FALSE]
    y <- y[complete_rows]
    if (nrow(A) < 10) next

    # Clustering
    founder_dist <- stats::dist(t(A))
    hc <- stats::hclust(founder_dist, method = "complete")
    groups <- stats::cutree(hc, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) cat(sprintf("Window %d: %d SNPs, %d groups\n", window_size, nrow(A), n_groups))

    # Require improvement over previous step
    if (final_result != NULL && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("  No improvement; continue\n")
      next
    }

    previous_n_groups <- n_groups

    # Constraints
    n_founders <- ncol(A)
    E <- matrix(1, nrow = 1, ncol = n_founders)  # sum-to-one
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("  Added %d accumulated constraints\n", nrow(accumulated_constraints)))
    }

    # Solve LSEI
    res <- tryCatch(
      limSolve::lsei(A = A, B = y, E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("  LSEI failed; continue\n")
      next
    }

    # Build constraints from current grouping
    current_constraints <- NULL
    current_constraint_values <- NULL
    for (cid in unique(groups)) {
      idx <- which(groups == cid)
      if (length(idx) > 1) {
        row <- rep(0, n_founders); row[idx] <- 1
        current_constraints <- rbind(current_constraints, row)
        current_constraint_values <- c(current_constraint_values, sum(res$X[idx]))
      } else {
        row <- rep(0, n_founders); row[idx] <- 1
        current_constraints <- rbind(current_constraints, row)
        current_constraint_values <- c(current_constraint_values, res$X[idx])
      }
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
    if (verbose >= 2) cat("  Stored result\n")
    if (n_groups == n_founders) {
      if (verbose >= 2) cat("  All founders distinguished; stop\n")
      break
    }
  }

  # Prepare outputs to mirror production
  if (!is.null(final_result)) {
    haps <- final_result$X
    names(haps) <- founder_names
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      err <- final_result$covar
    } else {
      err <- matrix(NA, length(founder_names), length(founder_names))
      rownames(err) <- founder_names
      colnames(err) <- founder_names
    }
  } else {
    haps <- rep(NA_real_, length(founder_names)); names(haps) <- founder_names
    err <- matrix(NA, length(founder_names), length(founder_names))
    rownames(err) <- founder_names; colnames(err) <- founder_names
  }

  list(Groups = groups, Haps = haps, Err = err, Names = founder_names)
}

# Centered-window variant that takes the full matrices once, then grows windows around pos
## Parameters:
## - founder_matrix: [num_snps × num_founders]
## - sample_freqs: [num_snps]
## - founder_names: length num_founders
## - positions: numeric vector [num_snps] genomic positions (bp)
## - pos: numeric, center position (bp)
## - window_sizes_bp: numeric vector of bp sizes, e.g., c(10000, 20000, ...)
## - h_cutoff, verbose: as above
## Returns: list(Groups, Haps, Err, Names)
estimate_haplotypes_centered <- function(
  founder_matrix,
  sample_freqs,
  founder_names,
  positions,
  pos,
  window_sizes_bp = c(10000, 20000, 50000, 100000, 200000, 500000),
  h_cutoff = 4,
  verbose = 0
) {
  stopifnot(is.matrix(founder_matrix))
  stopifnot(length(sample_freqs) == nrow(founder_matrix))
  stopifnot(length(founder_names) == ncol(founder_matrix))
  stopifnot(length(positions) == nrow(founder_matrix))

  final_result <- NULL
  final_n_groups <- 0
  previous_n_groups <- 0
  groups <- rep(1, length(founder_names))

  accumulated_constraints <- NULL
  accumulated_constraint_values <- NULL

  for (ws in window_sizes_bp) {
    left_bp  <- pos - ws/2
    right_bp <- pos + ws/2
    in_window <- positions >= left_bp & positions <= right_bp
    if (!any(in_window)) next

    A <- founder_matrix[in_window, , drop = FALSE]
    y <- sample_freqs[in_window]
    complete_rows <- stats::complete.cases(A) & !is.na(y)
    A <- A[complete_rows, , drop = FALSE]
    y <- y[complete_rows]
    if (nrow(A) < 10) next

    founder_dist <- stats::dist(t(A))
    hc <- stats::hclust(founder_dist, method = "complete")
    groups <- stats::cutree(hc, h = h_cutoff)
    n_groups <- length(unique(groups))
    if (verbose >= 2) cat(sprintf("Window %dbp: %d SNPs, %d groups\n", ws, nrow(A), n_groups))

    if (!is.null(final_result) && n_groups <= previous_n_groups) {
      if (verbose >= 2) cat("  No improvement; continue\n")
      next
    }
    previous_n_groups <- n_groups

    n_founders <- ncol(A)
    E <- matrix(1, nrow = 1, ncol = n_founders)
    F <- 1.0
    if (!is.null(accumulated_constraints)) {
      E <- rbind(E, accumulated_constraints)
      F <- c(F, accumulated_constraint_values)
      if (verbose >= 2) cat(sprintf("  Added %d accumulated constraints\n", nrow(accumulated_constraints)))
    }

    res <- tryCatch(
      limSolve::lsei(A = A, B = y, E = E, F = F, G = diag(n_founders), H = matrix(rep(0.0003, n_founders)), fulloutput = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || res$IsError != 0) {
      if (verbose >= 2) cat("  LSEI failed; continue\n")
      next
    }

    current_constraints <- NULL
    current_constraint_values <- NULL
    for (cid in unique(groups)) {
      idx <- which(groups == cid)
      if (length(idx) > 1) {
        row <- rep(0, n_founders); row[idx] <- 1
        current_constraints <- rbind(current_constraints, row)
        current_constraint_values <- c(current_constraint_values, sum(res$X[idx]))
      } else {
        row <- rep(0, n_founders); row[idx] <- 1
        current_constraints <- rbind(current_constraints, row)
        current_constraint_values <- c(current_constraint_values, res$X[idx])
      }
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
    if (verbose >= 2) cat("  Stored result\n")
    if (n_groups == n_founders) {
      if (verbose >= 2) cat("  All founders distinguished; stop\n")
      break
    }
  }

  if (!is.null(final_result)) {
    haps <- final_result$X
    names(haps) <- founder_names
    if ("covar" %in% names(final_result) && !is.null(final_result$covar)) {
      err <- final_result$covar
    } else {
      err <- matrix(NA, length(founder_names), length(founder_names))
      rownames(err) <- founder_names
      colnames(err) <- founder_names
    }
  } else {
    haps <- rep(NA_real_, length(founder_names)); names(haps) <- founder_names
    err <- matrix(NA, length(founder_names), length(founder_names))
    rownames(err) <- founder_names; colnames(err) <- founder_names
  }

  list(Groups = groups, Haps = haps, Err = err, Names = founder_names)
}

# Adapter to match production-like signature when desired
# Signature mirrors: estimate_haplotypes_list_format(pos, sample_name, df3, founders, ...)
# This adapter expects df3 to be in the tidy format with columns POS, name, freq where
# rows include all founders (names in `founders`) and the given sample.
# Window selection here is PREFIX-based over POS (assuming POS == 1..N). To switch back
# to center-based windows, replace the window extraction with a [center -/+ window] slice
# on POS and pass the resulting matrix/vector to estimate_haplotypes_prefix().
estimate_haplotypes_prefix_adapter <- function(
  pos, sample_name, df3, founders, h_cutoff,
  method = "adaptive", window_size_bp = NULL, chr = "chr2R", verbose = 0
) {
  # Build founder matrix and sample vector from df3 in POS order
  # Assumes POS is dense 1..N in simulated setting
  wide <- df3[, c("POS", "name", "freq")]
  wide <- stats::na.omit(wide)
  # founder matrix
  founder_wide <- reshape2::dcast(wide[wide$name %in% founders, ], POS ~ name, value.var = "freq")
  founder_wide <- founder_wide[order(founder_wide$POS), ]
  A <- as.matrix(founder_wide[, founders, drop = FALSE])
  # sample vector
  sample_wide <- reshape2::dcast(wide[wide$name == sample_name, ], POS ~ name, value.var = "freq")
  sample_wide <- sample_wide[order(sample_wide$POS), ]
  # align by POS
  common_pos <- intersect(founder_wide$POS, sample_wide$POS)
  founder_wide <- founder_wide[match(common_pos, founder_wide$POS), ]
  sample_wide <- sample_wide[match(common_pos, sample_wide$POS), ]
  A <- as.matrix(founder_wide[, founders, drop = FALSE])
  y <- as.numeric(sample_wide[[sample_name]])
  # Delegate to prefix estimator
  estimate_haplotypes_prefix(A, y, founders, h_cutoff = h_cutoff, verbose = verbose)
}

## No demo block: this module is called by the external wrapper that simulates data.


