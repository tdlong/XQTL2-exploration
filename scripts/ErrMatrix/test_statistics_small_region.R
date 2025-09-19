#!/usr/bin/env Rscript

# Simplified statistical test for small region
# Based on haps2scan functions but adapted for our test region

library(tidyverse)
library(limSolve)
library(abind)

# Load the reshaped data
reshaped_file <- "process/ZINC2_h10/adapt_h10/R.haps.chr3R.out.rds"
reshaped_data <- readRDS(reshaped_file)

# Define test region
test_positions <- seq(20000000, 20200000, by = 10000)

# Filter to test region
test_data <- reshaped_data %>%
  filter(pos %in% test_positions)

cat("=== STATISTICAL TESTING FOR SMALL REGION ===\n")
cat("Testing", nrow(test_data), "positions in region 20,000,000-20,200,000\n\n")

# Get sample names from the first position
sample_names <- test_data$sample[[1]]
cat("Sample names:", paste(sample_names, collapse = ", "), "\n\n")

# Create a simple design matrix for our samples
# Based on sample names: Rep01_W_M, Rep01_Z_M, etc.
design_df <- data.frame(
  bam = sample_names,
  TRT = ifelse(str_detect(sample_names, "_W_"), "W", "Z"),
  REP = str_extract(sample_names, "Rep\\d+"),
  REPrep = str_extract(sample_names, "Rep\\d+"),
  Num = 100,  # Assume 100 flies per sample
  Proportion = ifelse(str_detect(sample_names, "_Z_"), 0.5, NA)  # 50% selection for Z treatment
) %>%
  filter(str_detect(bam, "_M"))  # Only males

cat("Design matrix:\n")
print(design_df)

# Copy the statistical functions from scan_functions.R
average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  log_det <- sum(log(positive_eigenvalues))  
  n_positive <- length(positive_eigenvalues)  
  log_avg_var <- log_det / n_positive  
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

mn.covmat <- function(p, n, min.p = 0) {
  p[p < min.p] <- min.p
  p <- p / sum(p)
  mat <- -tcrossprod(p)
  diag(mat) <- p * (1 - p)
  mat <- mat / n
  mat
}

wald.test3 <- function(p1, p2, covar1, covar2, nrepl = 1, N1 = NA, N2 = NA) {
  if (nrepl > 1) {
    N1.eff <- rep(NA, nrepl)
    N2.eff <- rep(NA, nrepl)
    lp1 <- length(p1[1,])
    cv1 <- array(NA, c(lp1, lp1, nrepl))
    cv2 <- array(NA, c(lp1, lp1, nrepl))
    for (i in 1:nrepl) {
      covmat1 <- mn.covmat((N1[i] * p1[i,] + N2[i] * p2[i,]) / (N1[i] + N2[i]), 2 * N1[i])
      covmat2 <- mn.covmat((N1[i] * p1[i,] + N2[i] * p2[i,]) / (N1[i] + N2[i]), 2 * N2[i])
      N1.eff[i] <- sum(diag(covmat1)) * 4 * N1[i]^2 / (sum(diag(covmat1)) * 2 * N1[i] + 2 * N1[i] * sum(diag(covar1[,,i])))
      N2.eff[i] <- sum(diag(covmat2)) * 4 * N2[i]^2 / (sum(diag(covmat2)) * 2 * N2[i] + 2 * N2[i] * sum(diag(covar2[,,i])))
      cv1[,,i] <- (covmat1 + covar1[,,i]) * (N1.eff[i])^2
      cv2[,,i] <- (covmat2 + covar2[,,i]) * (N2.eff[i])^2
    }
    p1 <- N1.eff %*% p1 / sum(N1.eff)
    p2 <- N2.eff %*% p2 / sum(N2.eff)
    covar1 <- rowSums(cv1, dims = 2) / sum(N1.eff)^2
    covar2 <- rowSums(cv2, dims = 2) / sum(N2.eff)^2
  } else {
    covmat1 <- mn.covmat((N1 * p1 + N2 * p2) / (N1 + N2), 2 * N1)
    covmat2 <- mn.covmat((N1 * p1 + N2 * p2) / (N1 + N2), 2 * N2)
    covar1 <- covar1 + covmat1
    covar2 <- covar2 + covmat2
  }
  
  df <- length(p1) - 1
  covar <- covar1 + covar2
  eg <- eigen(covar)
  ev <- eg$vectors[,1:df]
  eval <- eg$values[1:df]
  trafo <- diag(1/sqrt(eval)) %*% t(ev)
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval <- exp(pchisq(tstat, df, lower.tail = FALSE, log.p = TRUE))
  list(wald.test = tstat, p.value = pval, avg.var = average_variance(covar)$avg_var)
}

pseudoN.test <- function(p1, p2, covar1, covar2, nrepl, N1, N2) {
  pseudoN_C <- rep(NA, nrepl)
  pseudoN_Z <- rep(NA, nrepl)
  for(i in 1:nrepl) {
    pseudoN_C[i] <- (2 * N1[i] * sum(p1[i,] * (1-p1[i]))) / (2 * N1[i] * sum(diag(covar1[,,i])) + sum(p1[i,] * (1-p1[i])))
    pseudoN_Z[i] <- (2 * N2[i] * sum(p2[i,] * (1-p2[i]))) / (2 * N2[i] * sum(diag(covar2[,,i])) + sum(p2[i,] * (1-p2[i])))
  }
  Count1 <- round(p1 * pseudoN_C, 0)
  Count2 <- round(p2 * pseudoN_Z, 0)
  lowCountFounder <- apply(rbind(Count1, Count2), 2, sum)
  if(sum(lowCountFounder >= 5) < 2) {
    log10p <- NA
  } else {
    Count1 <- Count1[,lowCountFounder >= 5]
    Count2 <- Count2[,lowCountFounder >= 5]
    if(nrepl == 1) {
      out <- chisq.test(rbind(Count1, Count2), correct = TRUE)
    } else {
      nF <- ncol(Count1)
      tdf <- data.frame(Count = c(as.numeric(t(Count1)), as.numeric(t(Count2))),
                       founder = rep(1:nF, 2*nrepl),
                       TRT = c(rep(1, nF*nrepl), rep(2, nF*nrepl)),
                       REP = c(rep(1:nrepl, each = nF), rep(1:nrepl, each = nF)))
      D.x <- xtabs(Count ~ founder + TRT + REP, data = tdf)
      out <- mantelhaen.test(D.x, correct = TRUE)
    }
    log10p <- -log10(out$p.value)
  }
  log10p
}

# Main testing function - Wald test only
doscan2 <- function(df, chr, Nfounders) {
  sexlink <- 1
  if(chr == "chrX") sexlink <- 0.75
  
  # Extract data for this position
  sample_names <- df$sample[[1]]
  groups <- df$Groups[[1]]
  haps <- df$Haps[[1]]
  err <- df$Err[[1]]
  names <- df$Names[[1]]
  
  # Create a data frame for this position
  df2 <- data.frame(
    sample = sample_names,
    Groups = groups,
    Haps = haps,
    Err = err,
    Names = names
  ) %>%
    left_join(design_df, join_by(sample == bam)) %>%
    filter(!is.na(TRT))
  
  # Check if all founders are discernable
  allFounders <- max(unlist(groups))
  
  ll <- list(Wald_log10p = NA, avg.var = NA)
  if(allFounders != Nfounders) return(ll)
  
  # Collapse replicates
  df3 <- df2 %>%
    select(-Groups) %>%
    group_by(TRT, REP) %>%
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`) / length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`) / length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink * mean(Num)) %>%
    rename(Haps = Haps_mean, Num = Num_mean, Err = Err_mean)
  
  # Extract data for testing
  p1 <- df3 %>% filter(TRT == "W") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p1) <- NULL
  p2 <- df3 %>% filter(TRT == "Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 <- do.call(abind, c(df3 %>% filter(TRT == "W") %>% pull(Err), along = 3))
  covar2 <- do.call(abind, c(df3 %>% filter(TRT == "Z") %>% pull(Err), along = 3))
  nrepl <- df3 %>% filter(TRT == "W") %>% nrow()
  N1 <- df3 %>% filter(TRT == "W") %>% pull(Num)
  N2 <- df3 %>% filter(TRT == "Z") %>% pull(Num)
  
  # Run Wald test only
  wt <- wald.test3(p1, p2, covar1, covar2, nrepl, N1, N2)
  Wald_log10p <- -log10(wt$p.value)
  
  ll <- list(Wald_log10p = Wald_log10p, avg.var = wt$avg.var)
  ll
}

# Run the tests
Nfounders <- length(test_data$Groups[[1]][[1]])

results <- test_data %>%
  group_by(CHROM, pos) %>%
  nest() %>%
  mutate(out = map2(data, CHROM, doscan2, Nfounders = Nfounders)) %>%
  unnest_wider(out)

# Display results
cat("=== WALD TEST RESULTS ===\n")
print(results %>% select(pos, Wald_log10p, avg.var), n = Inf)

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Positions tested:", nrow(results), "\n")
cat("Successful Wald tests:", sum(!is.na(results$Wald_log10p)), "\n")
cat("Mean Wald log10p:", round(mean(results$Wald_log10p, na.rm = TRUE), 2), "\n")
cat("Max Wald log10p:", round(max(results$Wald_log10p, na.rm = TRUE), 2), "\n")
cat("Min Wald log10p:", round(min(results$Wald_log10p, na.rm = TRUE), 2), "\n")
