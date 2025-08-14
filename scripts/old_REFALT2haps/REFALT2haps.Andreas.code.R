nrow_subset = function(spotsdf, df3){
        df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                nrow()
        }




est_hap = function(spotsdf, df3){
        # spotsdf = spots$data[[1]]  testing
        temp_mat = df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                select(-c(CHROM,N)) %>%
                pivot_wider(names_from=name, values_from=freq) %>%
                pivot_longer(!c("POS",matches(founders)),names_to = "sample", values_to = "freq") %>%
                select(-POS)

        sample_mat = temp_mat %>%
                group_by(sample) %>%
                nest() %>%
                mutate(haps=map(data,est_hap2)) %>%
                select(-data) %>%
                unnest_wider(haps)

        sample_mat
        }
        
# function for estimating the haplotype for a given sample, called from est_hap
# returns Groups from cuttree, hap freq estimates, errors by founder
est_hap2 = function(sampdf){

        # sampdf = sample_mat$data[[1]]  # testing
        founder_mat = sampdf %>% select(matches(founders))
        Y = sampdf$freq
        good = !is.na(Y)
        Y = Y[good]
        founder_mat = founder_mat[good,] 
        m_founder_mat = as.matrix(founder_mat)
        
        # Sanity checks for bad estimation space
        n_snps <- nrow(m_founder_mat)
        n_founders <- ncol(m_founder_mat)
        
        # Check 1: Too few SNPs relative to founders (rule of thumb: need at least 3x)
        if (n_snps < n_founders * 3) {
                warning("Bad estimation space: ", n_snps, " SNPs for ", n_founders, " founders (need at least ", n_founders * 3, ")")
                return(list(Groups = rep(1, n_founders), 
                           Haps = rep(NA, n_founders), 
                           Err = matrix(NA, n_founders, n_founders), 
                           Names = names(founder_mat)))
        }
        
        # Check 2: Matrix condition number (numerical stability)
        if (n_snps >= n_founders) {
                condition_num <- kappa(m_founder_mat)
                if (condition_num > 1e10) {
                        warning("Bad estimation space: Matrix condition number too high (", format(condition_num, scientific = TRUE), ")")
                        return(list(Groups = rep(1, n_founders), 
                                   Haps = rep(NA, n_founders), 
                                   Err = matrix(NA, n_founders, n_founders), 
                                   Names = names(founder_mat)))
                }
        }
        
        # Check 3: Effective rank (how many founders are actually distinguishable)
        if (n_snps >= n_founders) {
                svd_result <- svd(m_founder_mat)
                effective_rank <- sum(svd_result$d > 1e-6)
                if (effective_rank < n_founders * 0.7) {
                        warning("Bad estimation space: Effective rank too low (", effective_rank, " for ", n_founders, " founders)")
                        return(list(Groups = rep(1, n_founders), 
                                   Haps = rep(NA, n_founders), 
                                   Err = matrix(NA, n_founders, n_founders), 
                                   Names = names(founder_mat)))
                }
        }
        
        Groups = cuttree(hclust(dist(t(m_founder_mat))),h=h_cutoff)
        d = ncol(m_founder_mat)         
        out = lsei(A=m_founder_mat,B=Y, E=t(matrix(rep(1,d))),F=1,G=diag(d),H=matrix(rep(0.0003,d)),verbose=TRUE,fulloutput=TRUE)
        Haps = out$X
        Err = out$cov
        list(Groups=Groups,Haps=Haps,Err=Err,Names=names(Haps))
        }

df = lazy_dt(read.table(filein,header=TRUE))
df2 = df %>%
	pivot_longer(c(-CHROM,-POS), names_to = "lab", values_to = "count") %>%
	mutate(RefAlt = str_sub(lab,1,3)) %>%
	mutate(name = str_sub(lab,5)) %>%
	select(-lab) %>%
#	separate(lab, c("RefAlt", "name"), "_", extra = "merge") %>%
	pivot_wider(names_from = RefAlt, values_from = count) %>%
	mutate(freq = REF/(REF+ALT), N = REF+ALT) %>%
	select(-c("REF","ALT")) %>%
	as_tibble()

rm(df)
cat("df2 is now made\n")

# identify SNPs that are NOT problematic in the set of founders
good_SNPs = df2 %>%
	filter(name %in% founders) %>%
	group_by(CHROM,POS) %>%
	summarize(zeros=sum(N==0),notfixed=sum(N!=0 & freq > 0.03 & freq < 0.97),informative=(sum(freq)>0.05 | sum(freq) < 0.95)) %>%
	ungroup() %>%
	filter(zeros==0 & notfixed==0 & informative=="TRUE") %>%
	select(c(CHROM,POS))

# now subset the entire dataset for the good SNPs only
df3 = good_SNPs %>%
	left_join(df2, multiple = "all")

rm(df2)
cat("df3 is now made\n")
saveRDS(df3, file = rdsfile)

# df3 = readRDS(rdsfile)
# spots are the locations at which we will estimate haplotypes
# every <step> bp (i.e., 10kb) on the step (i.e, 0, 10, 20, ... kb)
# I define a window +/- size on those steps, and fix the ends
minpos = min(df3$POS)
maxpos = max(df3$POS)
myseq = seq(0,maxpos,step)
myseq = myseq[myseq > minpos + size & myseq < maxpos - size]
spots = data.frame(CHROM=rep(mychr,length(myseq)), pos=myseq, start=myseq-size, end=myseq+size)
# get rid of windows with fewer than 50 SNPs
# i.e., <50 SNPs in 100kb is pretty strange
UU = unique(df3$POS)
spots = spots %>%
                rowwise() %>%
                mutate(NN = sum(start < UU) - sum(end < UU)) %>%
                filter(NN >= 50) %>%
                select(-NN)

# this is the actual scan
# I guess this could also be slow...
# as it runs for L loci X S samples
spots2 = spots %>%
        group_nest(row_number()) %>%
        mutate(out = map(data,est_hap,df3)) %>%
		unnest(data) %>%
        select(-c(start,end)) %>%
        select(-`row_number()`) %>%
        unnest_wider(out) 

saveRDS(spots2,file = fileout)

