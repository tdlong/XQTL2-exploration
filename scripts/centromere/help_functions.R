library(tidyverse)
library(abind)

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

	## these summaries of the data are pretty useful for tests
Wald_wrapper = function(df3){
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

	wt=wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2)
	Wald_log10p = -log10(wt$p.value)
    return(Wald_log10p)
}

# Centromere-specific wrapper that works with the centromere data structure
Wald_wrapper_centromere = function(df3){
	# Extract haplotype frequencies for each treatment group
	p1 = df3 %>% filter(TRT=="W") %>% pull(Haps) %>% map_dfr(~as_tibble(t(.x))) %>% as.matrix()
	p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% map_dfr(~as_tibble(t(.x))) %>% as.matrix()
	
	# Extract error covariance matrices
	covar1 = df3 %>% filter(TRT=="W") %>% pull(Err) %>% abind::abind(along = 3)
	covar2 = df3 %>% filter(TRT=="Z") %>% pull(Err) %>% abind::abind(along = 3)
	
	# Sample sizes (assuming equal for now - you may need to adjust this)
	nrepl = df3 %>% filter(TRT=="W") %>% nrow()
	nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
	
	# For centromere data, we don't have Num column, so use a default or calculate from data
	# You may need to adjust this based on your actual sample sizes
	N1 = rep(1, nrepl)  # Default sample size - adjust as needed
	N2 = rep(1, nrepl)  # Default sample size - adjust as needed

	wt = wald.test3(p1, p2, covar1, covar2, nrepl, N1, N2)
	Wald_log10p = -log10(wt$p.value)
    return(Wald_log10p)
}

# Calculate average haplotype frequency difference between treatments
calculate_haplotype_difference <- function(data) {
  # Calculate average haplotype frequencies for each treatment
  avg_C <- data %>% 
    filter(TRT == "C") %>% 
    pull(Haps) %>% 
    map_dfr(~as_tibble(t(.x))) %>% 
    summarise_all(mean) %>% 
    as.numeric()
  
  avg_Z <- data %>% 
    filter(TRT == "Z") %>% 
    pull(Haps) %>% 
    map_dfr(~as_tibble(t(.x))) %>% 
    summarise_all(mean) %>% 
    as.numeric()
  
  # Calculate difference (Z - C)
  Dhap <- avg_Z - avg_C
  return(Dhap)
}

# Complete centromere analysis function
run_centromere_analysis <- function(results_file = "combined_centromere_all_results.RDS", 
                                   info_file = "info.ZINC2.txt",
                                   output_file = "Centromere_Wald_testing.RDS") {
  
  cat("=== CENTROMERE WALD TESTING ANALYSIS ===\n")
  cat("Loading data...\n")
  
  # Read the info file
  info_data <- read_tsv(info_file)
  cat("✓ Loaded info file:", info_file, "\n")
  
  # Run the complete analysis
  cat("Running analysis...\n")
  out <- readRDS(results_file) %>%
    separate(sample, into = c("rep", "TRT", "sex"), sep = "_", remove = FALSE) %>%
    left_join(info_data, by = c("sample" = "bam")) %>%
    mutate(TRT = ifelse(TRT == "W", "C", TRT)) %>%
    group_by(CHROM, pos, sex) %>%
    nest() %>%
    mutate(
      wald_results = map(data, Wald_wrapper),
      Dhap = map(data, calculate_haplotype_difference)
    )
  
  cat("✓ Analysis complete\n")
  cat("Results summary:\n")
  print(out)
  
  # Save results
  saveRDS(out, output_file)
  cat("✓ Results saved to:", output_file, "\n")
  
  return(out)
}

# Display Dhap results as percentage table
display_dhap_table <- function(results, founder_names = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "AB8")) {
  # Extract Dhap data and convert to percentages
  dhap_table <- results %>%
    select(CHROM, pos, sex, Dhap) %>%
    mutate(
      Dhap_pct = map(Dhap, ~ round(.x * 100, 2))
    ) %>%
    select(-Dhap) %>%
    unnest(Dhap_pct) %>%
    mutate(founder = rep(founder_names, nrow(results))) %>%
    pivot_wider(names_from = founder, values_from = Dhap_pct) %>%
    select(CHROM, pos, sex, all_of(founder_names))
  
  return(dhap_table)
}

