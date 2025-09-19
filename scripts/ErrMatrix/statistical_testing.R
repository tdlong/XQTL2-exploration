#!/usr/bin/env Rscript

# Statistical testing script for haplotype data
# Usage: Rscript statistical_testing.R <data_dir> <chromosome> <start_pos> <end_pos> <design_file>
# Example: Rscript statistical_testing.R process/ZINC2_h10/adapt_h10 chr3R 20000000 20200000 /dfs7/adl/tdlong/fly_pool/XQTL2/helpfiles/ZINC2/Zinc2.test.M.txt

library(tidyverse)
library(limSolve)
library(abind)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript statistical_testing.R <data_dir> <chromosome> <start_pos> <end_pos> <design_file>")
}

mydir <- args[1]
mychr <- args[2]
start_pos <- as.numeric(args[3])
end_pos <- as.numeric(args[4])
design_file <- args[5]

# Load design file and data
design.df <- read.table(design_file)
source("scripts/haps2scan/scan_functions.R")
filein <- paste0(mydir, "/R.haps.", mychr, ".out.rds")

cat("=== STATISTICAL TESTING PARAMETERS ===\n")
cat("Data directory:", mydir, "\n")
cat("Chromosome:", mychr, "\n")
cat("Position range:", start_pos, "-", end_pos, "\n")
cat("Design file:", design_file, "\n")
cat("Input file:", filein, "\n\n")

library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)


average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  # Calculate eigenvalues
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  # Filter out eigenvalues that are effectively zero or negative
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  # Calculate the product of positive eigenvalues
  log_det <- sum(log(positive_eigenvalues))  
  # Use the number of positive eigenvalues for the root
  n_positive <- length(positive_eigenvalues)  
  # Calculate log of average variance
  log_avg_var <- log_det / n_positive  
  # Convert back to original scale
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

wald.test3 = function(p1,p2,covar1,covar2,nrepl=1,N1=NA,N2=NA){
    
    # Wald test for multinomial frequencies
    # if nrepl = 1: (one replicate, analogous to chi square):
    #  p1 and p2 are vectors of relative frequencies to be compared
    # covar1 and covar2 are the reconstruction error 
    # covariance matrices from limSolve
    # the sampling covariance matrices are generated within limSolve
    # if nrepl > 1 (multiple replicates, analogous to CMH):
    #   p1 and p2 are matrices, each row is frequency vector for one replicate
    # covar1 and covar2 are tensors (3-dimensional arrays, third dimension 
    #  denotes replicate) for the linSolve covariance matrices
    # N1 (initial) and N2 (after treatment) 
    # are sample sizes, they are vectors when there is more than one replicate
    # N1[i], N2[i] are then for replicate i
    if (nrepl>1){
      N1.eff=rep(NA,nrepl)
      N2.eff=rep(NA,nrepl)
      lp1 = length(p1[1,])
      cv1=array(NA,c(lp1,lp1,nrepl))
      cv2=array(NA,c(lp1,lp1,nrepl))
      for (i in 1:nrepl){
          
        covmat1  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N1[i])
        covmat2  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N2[i])
    
        N1.eff[i] = sum(diag(covmat1))*4*N1[i]^2/(sum(diag(covmat1))*2*N1[i]+2*N1[i]*sum(diag(covar1[,,i])) )
        N2.eff[i] = sum(diag(covmat2))*4*N2[i]^2/(sum(diag(covmat2))*2*N2[i]+2*N2[i]*sum(diag(covar2[,,i])) )
        cv1[,,i]= (covmat1 + covar1[ , ,i])  * (N1.eff[i])^2
        cv2[,,i]= (covmat2 + covar2[ , ,i])  * (N2.eff[i])^2
        
      }
        
      p1 = N1.eff %*% p1 / sum(N1.eff)
      p2 = N2.eff %*% p2 / sum(N2.eff)
     
      covar1= rowSums(cv1, dims = 2) / sum(N1.eff)^2
      covar2= rowSums(cv2, dims = 2) / sum(N2.eff)^2
     # browser()
    }
    else {
      covmat1  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N1)
      covmat2  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N2)
      covar1 = covar1 + covmat1
      covar2 = covar2 + covmat2
    }
  
  df = length(p1)-1
  covar=covar1+covar2
  eg <- eigen(covar)
  # remove last eigenvector which corresponds to eigenvalue zero
  ev <- eg$vectors[,1:df]
  eval <- eg$values[1:df]
  trafo<-diag(1/sqrt(eval)) %*% t(ev) 
  # set extremely small values to zero
  #new.covar[new.covar < 10^-9]=0
  p1= as.vector(p1); p2=as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval<- exp(pchisq(tstat,df,lower.tail=FALSE,log.p=TRUE))
  list(wald.test=tstat, p.value=pval, avg.var=average_variance(covar)$avg_var)
}

mn.covmat= function(p,n,min.p=0){
  # generate multinomial covariance matrix
  # p is vector of multinomial relative frequencies
  # n is sample size
  # compute covariance matrix for relative frequencies, for absolute frequencies multiply by n^2
  # if min.p >0, then values of p smaller than min.p are set to min.p and the resulting vector is rescaled.
  p[p<min.p] = min.p; p=p/sum(p)
  mat = - tcrossprod(p)
  diag(mat) = p*(1-p)
  mat = mat/n
  mat
}
doscan2 = function(df,chr,Nfounders){
	sexlink = 1
	if(chr=="chrX"){ sexlink=0.75 }

	# I tested with xx2$data[[1]]
	df2 = df %>%
		unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
		left_join(design.df, join_by(sample==bam)) %>%
		filter(!is.na(TRT))
	
	# only analyze data for which all founders are discernable..
	allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	

	ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
	if(allFounders!=Nfounders){ return(ll) }

	##  now cases where all founders are OK
	##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
	##  is value in keeping dataframe columns as lists...
	df3 = df2 %>%
		select(-Groups) %>%
		group_by(TRT,REP) %>%	
		summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
			Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
			Names = list(first(Names)),
			Num_mean = sexlink*mean(Num)) %>%
		rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)

	## these summaries of the data are pretty useful for tests
	p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
	row.names(p1) <- NULL
	p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
	row.names(p2) <- NULL
	covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
	covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
	nrepl = df3 %>% filter(TRT=="C") %>% nrow()
	nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
	N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
	N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)

	# Store haplotype frequencies for analysis
	hap_freqs_C = as.vector(p1)
	hap_freqs_Z = as.vector(p2)
	hap_diff = hap_freqs_C - hap_freqs_Z
	
	# Store error variances for analysis
	err_var_C = diag(covar1[,,1])
	err_var_Z = diag(covar2[,,1])
	err_diff = err_var_C - err_var_Z

	wt=wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2)
	Wald_log10p = -log10(wt$p.value)
#	Pseu_log10p = pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2)

#	af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
#	temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
#	Falc_H2 = temp$Falconer_H2
#	Cutl_H2 = temp$Cutler_H2

	ll = list(Wald_log10p = Wald_log10p, 
	          hap_freqs_C = list(hap_freqs_C),
	          hap_freqs_Z = list(hap_freqs_Z), 
	          hap_diff = list(hap_diff),
	          err_var_C = list(err_var_C),
	          err_var_Z = list(err_var_Z),
	          err_diff = list(err_diff))
	ll
	}

xx1 = readRDS(filein)
Nfounders=length(xx1$Groups[[1]][[1]])
ProportionSelect = design.df %>% filter(TRT=="Z") %>% select(REP,Proportion) %>% arrange(REP)

bb1 = xx1 %>%
#	head(n=100) %>%
	filter(pos >= start_pos & pos <= end_pos) %>%  # User-specified region
	group_by(CHROM,pos) %>%
	nest() %>%
	mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
	unnest_wider(out)
bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)
bb3 = add_genetic(bb2)

# Display results for our test region
cat("=== WALD TEST RESULTS FOR TEST REGION ===\n")
print(bb2 %>% select(pos, Wald_log10p), n = Inf)

# Analyze haplotype frequency changes
cat("\n=== HAPLOTYPE FREQUENCY ANALYSIS ===\n")

# Extract haplotype frequencies and calculate changes
hap_analysis <- bb2 %>%
  filter(!is.na(Wald_log10p)) %>%
  arrange(pos) %>%
  select(pos, hap_freqs_C, hap_freqs_Z, hap_diff)

# Calculate changes between adjacent positions
hap_changes <- hap_analysis %>%
  mutate(
    hap_diff_prev = lag(hap_diff),
    hap_change = map2(hap_diff, hap_diff_prev, ~ abs(as.numeric(.x) - as.numeric(.y)))
  ) %>%
  filter(!is.na(hap_diff_prev)) %>%
  select(pos, hap_change)

cat("Haplotype treatment differences (C - Z) by position:\n")
for(i in 1:nrow(hap_analysis)) {
  cat("Position", hap_analysis$pos[i], ":", 
      paste(round(as.numeric(hap_analysis$hap_diff[[i]]), 4), collapse = " "), "\n")
}

cat("\nChanges in haplotype differences between adjacent positions:\n")
for(i in 1:nrow(hap_changes)) {
  cat("Position", hap_changes$pos[i], ":", 
      paste(round(as.numeric(hap_changes$hap_change[[i]]), 4), collapse = " "), "\n")
}

# Calculate summary statistics
all_hap_changes <- unlist(hap_changes$hap_change)
cat("\nHaplotype change summary:\n")
cat("Mean change:", round(mean(all_hap_changes), 4), "\n")
cat("SD change:", round(sd(all_hap_changes), 4), "\n")
cat("Max change:", round(max(all_hap_changes), 4), "\n")

cat("\n=== SUMMARY ===\n")
cat("Positions tested:", nrow(bb2), "\n")
cat("Successful Wald tests:", sum(!is.na(bb2$Wald_log10p)), "\n")
cat("Mean Wald log10p:", round(mean(bb2$Wald_log10p, na.rm = TRUE), 2), "\n")
cat("Max Wald log10p:", round(max(bb2$Wald_log10p, na.rm = TRUE), 2), "\n")
cat("Min Wald log10p:", round(min(bb2$Wald_log10p, na.rm = TRUE), 2), "\n")

