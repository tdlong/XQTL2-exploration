# XQTL2 Pipeline Scripts

This repository contains the bioinformatics pipeline scripts for XQTL (Experimental Quantitative Trait Locus) analysis in Drosophila.

**Note:** The XQTL analysis and visualization functions have been moved to the [XQTL2.Xplore R package](https://github.com/tdlong/XQTL2.Xplore), which provides comprehensive tools for data analysis and plotting.

## Script Organization

The scripts are organized into logical subfolders reflecting the pipeline flow:

### **Main Scripts** (`scripts/`)
- **Main Pipeline**: `haplotype_testing_from_table.sh` - Comprehensive wrapper for haplotype estimation and SNP imputation
- **Haplotype Estimation**: `REFALT2haps.FixedWindow.Single.R`, `REFALT2haps.AdaptWindow.Single.R` - Single-parameter estimators
- **SNP Imputation**: `euchromatic_SNP_imputation_single.R` - SNP frequency interpolation with founder order fix
- **Method Evaluation**: `evaluate_haplotype_methods.R` - Comprehensive comparison of different haplotype estimation methods

### **Pipeline Subfolders**
- **`raw2bam2REFALT/`** - Early pipeline: Raw data → BAM → REF/ALT counts
  - `fq2bam.sh`, `fq2bam.Oct28.sh`, `bam2bcf2REFALT.sh`
- **`haps2scan/`** - GWAS scanning and analysis
  - All `haps2scan.*` scripts, `scan_functions.R`, `XQTL_plotting_functions.R`, `concat_Chromosome_Scans.*`
- **`old_REFALT2haps/`** - Legacy Andreas scripts
  - `REFALT2haps.Andreas.*` scripts
- **`Heterozygosity_tests/`** - Heterozygosity analysis
  - `plot_heterozygosity.R`, `process_heterozygosity.R`, `bam2bcf2REFALT_het.sh`, `REFALT2HET.awk`
- **`debug_and_testing/`** - Debugging, testing, and diagnostic scripts (16 files)

## Current XQTL pipeline

## Get sequences from core...

```bash
# save the email from Yuzo as blah.txt and then extrat the files
# you will have to manually edit the resulting file a little
cat blah.txt | grep http | cut -f1 -d' ' | awk '{printf("wget %s\n",$0)}' >get_data.sh

# get_data.sh lines should look like this
wget https://hts.igb.uci.edu/tdlong24102845/xR019-L8-G2-P045-TGGCTATG-TTGTCAGC-READ1-Sequences.txt.gz
# make a directory to store your raw data
mkdir data/raw/Oct28_24
# cd to folder, put get_data.sh script in folder, then run it as a slurm job
sbatch get_data.sh
# you will have to add the 4 lines below to the top of the get_data.sh file
#  so it is a slurm script, for example something like the below would allow 
#  it to run, all my scripts are designed to run on a slurm cluster and have
#  lines like the below in them

#!/bin/bash
#SBATCH --job-name=getdata
#SBATCH -A ???        ## account to charge 
#SBATCH -p standard   ## partition/queue name
```
## Create a file with read name mapping to fastq files.  And map the reads to the reference genome

```bash
# first you need indexed genomes to align to
# my script expect them to be in ref/ and have the root filename dm6.fa
# you have to be careful to map your reads to the same genome as the founders!
# you can get an indexed reference genome from here
cp /dfs7/adl/tdlong/fly_pool/newpipeline_aging/XQTL_pipeline/ref/* ref/.
# here is what is now in your directory (standard indexing)
# the files are big enough to not be on git
ls ref
dm6.dict
dm6.fa
dm6.fa.amb
dm6.fa.ann
dm6.fa.bwt
dm6.fa.fai
dm6.fa.pac
dm6.fa.sa
dm6.fa.sizes
README_ref.txt
README.txt

# now you need a file that lets your name the bam files from the raw read barcodes
# There are many ways to do this, I map barcodes to read names
# In the file below the first column is F barcode, 2nd column the R barcode,
#  and 3rd is sample name (tab delimited)
# put this file in the helperfiles directory
# note the structure of the raw file names above, depending on how your
# raw files are named my scripts may need to be modified
cat helperfiles/readname.mapping.JUICE.txt
TGGCTATG	TTGTCAGC	JUICE1
GTCCTAGA	TTGTCAGC	JUICE2
ACTTGCCA	TTGTCAGC	JUICE3
TCTTCGTG	TTGTCAGC	JUICE4
TCCACTCA	TTGTCAGC	JUICE5
CTGTGCTT	TTGTCAGC	JUICE6

# now run the alignments (fq -> bam)
# you may need to define the dir where the raw data is and where you want
#  to write bams to in the script (data/raw/JUICE below)
# you will need to make the array job as big as the number of samples being aligned
# (the number of line in readname mapping file above)
mkdir data/bam/JUICE
NN=`wc -l helperfiles/readname.mapping.JUICE.txt | cut -f1 -d' '`
sbatch --array=1-$NN scripts/raw2bam2REFALT/fq2bam.sh helperfiles/readname.mapping.JUICE.txt data/raw/JUICE data/bam/JUICE 

# after it finishes (it could take overnight)
ls -alh data/bam/JUICE

14G  ... data/.../JUICE1.bam
18G  ... data/.../JUICE2.bam
16G  ... data/.../JUICE3.bam
14G  ... data/.../JUICE4.bam
9.0G ... data/.../JUICE5.bam
17G  ... data/.../JUICE6.bam
```
It may be helpful to look at the alignment script above. We align (bwa mem), sort, add readgroups, index.  The readgroups are important.  The file sizes of the bams are also important (as files below about 1G are likely of insufficient coverage or could indicate a failed library prep). 

## Go from bams to bcf to REFALT 
The scripts in this section produce files that tabulates counts of REF and ALT alleles (for well behaved SNPs) for each SNP and sample

During this step we create a file ("helpfiles/Oct28_24.bams") that provides paths to all the bam files.  This file is important as it contains paths to all your poolseq samples, plus the pre-aligned founder bams.  Since the pre-aligned bams are big, I just give paths to where I store them on hpc3, clearly this requires you are part of my "group" and can see them.  If you not doing this at UCI, then you would need to download these bam files from somewhere and they would take several TBs of space.  Note how I exclude small bams, ideally I reprep and resequence these samples to have a complete dataset.
```bash
# This could take 20+ hours
# I will put processed stuff here
mkdir process
mkdir process/JUICE
# add bam files, then add founders
find data/bam/JUICE -name "*.bam" -size +1G > helpfiles/JUICE.bams
cat helpfiles/founder.bams.txt | grep "B" >>helpfiles/JUICE.bams
# now generate the REFALT files
sbatch scripts/raw2bam2REFALT/bam2bcf2REFALT.sh helpfiles/JUICE.bams process/JUICE
```

## Edit JUICE_haplotype_parameters.R to reflect your data
This file contains a bunch of information that lets scripts that call haplotypes and run the "GWAS" scan without intervention.  So you will probably spend some time tweaking this file.  If you do sub-analyses (on subsets of your data) your would have multiple version of this file that are passed to scripts. The file is mostly comments!
```bash
##########	
#  R project specific parameters
#  haplotype.parameters.R
#  note that the number are here for when I move to 6 replicates
##########
# RGs = the readgroups add to the bams earlier 
# In my pipeline the RGs and bam file names (less the extension) are the same

# list of founders for the specific population of the experiment
# take from the RGs of the founder bams, and generally abbreviate the standard DSPR way
# different or custom founder populations would have different founders
# the software estimates the frequency of each founder haplotype in each pooled sample

founders=c("B1","B2","B3","B4","B5","B6","B7","AB8")

# list of samples to process
# used in REFALT2hap 
# the list of sample names to consider, they should match up with the 
# names given to the earlier fq2bam script, as these are taken from the RGs of the resulting bams
# the command line below will generate the same list the one's with RAFALT data at an earlier step
# mybams="data/bam/JUICE"
# mysize="1G"
# echo -n "names_in_bam=c(" && find $mybams -name "*.bam" -size +$mysize -print0 |\
# 	xargs -0 -n1 basename |\
#	sed 's/.*/"&"/' |\
#	tr '\n' ',' |\
#	sed 's/,$//' && echo ")"

names_in_bam=c("R1con","R1age","R2con","R2age","R3con","R3age","R4con","R4age","R5con","R5age","R6con","R6age")

# step_size (bp)
# haplotypes will be imputed every step/1000 kb at the kb (each @ 10,20,30,etc)

step = 10000

# +/- window_size (bp)
# the windows are centered on the steps above, but of width +/- size/1000 kb
# if Numflies is generally large (>>200) and seq coverage high I think 50kb is good here
# but this is for sure a tuning parameter in that bigger windows lead to better
# haplotype inference, but poorer localization.  Smaller windows lead to better
# localization (in theory), but that is offset by poor haplotype inference.  One
# indication the window is too small is noisy neighbouring -log10p values from
# the scan.

size = 50000

# tree cutoff height to claim founders cannot be distinguished from one another
# this parameter is somewhat mysterious if I am being honest
# its units are Euclidean distance = sqrt(sum((Fi-Fj)^2) between two founders
# so over a window of 500 SNPs a distance of <2.5 implies only
# 6 SNPs being different, or 1.2% of SNPs, pretty similar haplotypes
# a cutoff of 5 implies 25 SNPs different, or 5% of SNPs
# a problem is the distance is not corrected for the number of SNPs, so if a window
# has 1000 as opposed to 500 SNPs ... then a given distance implies fewer fixed
# difference SNPs. I would tend to leave this at 2.5 or perhaps 5.

h_cutoff=2.5

```
The bash command below can be useful for editting the haplotype parameters file, as in many cases only the list of bam files is changing (and perhaps the A vs B founders).  You only need to edit the path to the bams and perhaps the size cutoff for considering the bam file not a "redo" to get a list of bams.
```bash
mybams="data/bam/STARVE"
mysize="1G"
echo -n "names_in_bam=c(" && find $mybams -name "*.bam" -size +$mysize -print0 |\
 	xargs -0 -n1 basename |\
	sed 's/.bam//' |\
	sort |\
	sed 's/.*/"&"/' |\
	tr '\n' ',' |\
	sed 's/,$//' && echo ")"
```

## Haplotype Estimation and SNP Imputation (30-60min)

The main pipeline now uses a streamlined approach with a comprehensive wrapper that can test multiple parameter combinations and optionally run SNP imputation.

### Quick Start: Test chr2R with Multiple Parameters

Create a parameter table for testing different approaches:

```bash
# Create parameter table for chr2R testing
cat > helpfiles/haplotype_params.2R.tsv << 'EOF'
chr2R	fixed	20
chr2R	fixed	50
chr2R	fixed	100
chr2R	fixed	200
chr2R	fixed	500
chr2R	adaptive	4
chr2R	adaptive	6
chr2R	adaptive	8
chr2R	adaptive	10
EOF
```

Run the main pipeline (9 array jobs):

```bash
# Haplotype estimation only
sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE

# Haplotype estimation + SNP imputation
sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE yes
```

### What This Does

**Haplotype Estimation:**
- **Fixed window analysis** - Tests 20kb, 50kb, 100kb, 200kb, 500kb windows
- **Adaptive window analysis** - Tests h_cutoff values of 4, 6, 8, 10
- **Euchromatin filtering** - Automatically filters to euchromatic regions only
- **Parallel processing** - Each parameter combination runs as a separate array job

**Optional SNP Imputation:**
- **Interpolation** - Uses haplotype estimates to interpolate SNP frequencies
- **Founder order fix** - Includes the critical founder order bug fix
- **Euchromatin focus** - Only imputes SNPs in euchromatic regions

### Output Files

**Haplotype Results:**
- `fixed_window_*kb_results_chr2R.RDS` - Fixed window haplotype estimates
- `adaptive_window_h*_results_chr2R.RDS` - Adaptive window haplotype estimates

**SNP Imputation Results (if enabled):**
- `snp_imputation_fixed_*kb_chr2R.RDS` - Imputed SNP frequencies for fixed windows
- `snp_imputation_adaptive_h*_chr2R.RDS` - Imputed SNP frequencies for adaptive windows

### Method Evaluation

After running multiple parameter combinations, evaluate which method performs best:

```bash
# Evaluate all methods on chr2R
Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE 1000
```

**What the evaluation does:**
- **MSE Analysis** - Calculates Mean Squared Error between observed and imputed SNP frequencies
- **Coverage Metrics** - Measures proportion of positions with valid estimates
- **Regional Performance** - Analyzes MSE across chromosomal windows (default 1000kb)
- **Method Comparison** - Generates summary statistics and recommendations
- **Visualizations** - Creates comprehensive plots comparing all methods

**Evaluation Output:**
- `haplotype_evaluation_summary_chr2R.RDS` - Overall performance summary
- `haplotype_evaluation_detailed_chr2R.RDS` - Sample-by-sample metrics
- `haplotype_evaluation_regional_chr2R.RDS` - Regional performance analysis
- `haplotype_evaluation_plots_chr2R.png` - Comprehensive comparison plots

**Key Metrics:**
- **MSE/RMSE** - Accuracy of imputed frequencies
- **Coverage** - Proportion of positions with estimates
- **Correlation** - Correlation with observed frequencies
- **Regional Analysis** - Performance across chromosomal regions

### Euchromatin Boundaries

The pipeline automatically uses the correct euchromatin boundaries for each chromosome:
- **chr2L**: 82,455 - 22,011,009 bp
- **chr2R**: 5,398,184 - 24,684,540 bp  
- **chr3L**: 158,639 - 22,962,476 bp
- **chr3R**: 4,552,934 - 31,845,060 bp
- **chrX**: 277,911 - 22,628,490 bp

### Legacy Approach (Still Available)

For legacy compatibility, the old Andreas scripts are still available:

```bash
# Legacy haplotype estimation
sbatch scripts/old_REFALT2haps/REFALT2haps.Andreas.sh helpfiles/haplotype_parameters.R "process/Oct28_24"
```

## Run the scan (15min) -- the old way
```bash
# point to same folder as above for in files
# note libraries needed
# the last argument defines a folder for output (it is created by the script) and is inside the input folder
# There is a testing parameters file that can be changed for different analyses
sbatch scripts/haps2scan/haps2scan.Andreas.sh helperfiles/testing.parameters.A.R "process/Oct28_24" "TEST_A"
sbatch scripts/haps2scan/haps2scan.Andreas.sh helperfiles/testing.parameters.B.R "process/Oct28_24" "TEST_B"
```
Here is an example of the testing folder, these are painful to make..

```bash
##########	
#  R project specific parameters
#  testing.parameters.R
#  note that the number are here for when I move to 6 replicates
##########

# list of samples
# used in REFALT2hap 
# the list of sample names to consider, they should match up with the 
# names given to the earlier fq2bam, as these are taken from the RGs of the bam 
names_in_bam=c("R1con","R1age","R2con","R2age","R3con","R3age","R4con","R4age","R5con","R5age","R6con","R6age")

# note the naming convention has three fields
#  Con vs Treatment -- must have two levels
#  Replicate
#  Possibly replicate within replicate, often "1"
samples=c("Con_1_1","Age_1_1","Con_2_1","Age_2_1","Con_3_1","Age_3_1","Con_4_1","Age_4_1","Con_5_1","Age_5_1","Con_6_1","Age_6_1")

# Numflies
# The number of flies in each pool
Numflies = data.frame(pool=samples,Num=c(570,1177,520,814,610,482,580,997,640,542,610,647))

# Proportion of Flies selected per replicate
ProportionSelect = data.frame(REP=c(1,2,3,4,5,6),Proportion=c(0.113,0.087,0.040,0.080,0.045,0.053))

# Mapping of Treatments to Control versus Selected
# Prefixes must be mapped to C for controls or Z for selected
TreatmentMapping = data.frame(longTRT=c("Con","Age"),TRT=c("C","Z"))

```
## Run the scan (15min) -- the new way

The new way seems a little easier.  Like the old way there is a path to the input data and a new folder you define where the output goes (inside the input folder).  But with the new way, instead of a parameter file, you just point to a file that can be read into R via "read.table".  It is best to point to an object saved via "write.table" in R with no other switches.  

```bash
sbatch scripts/haps2scan/haps2scan.Apr2025.sh helpfiles/Oct28_24.testA.txt "process/Oct28_24" "TEST_A"

```
The R dataframe requires certain columns -- their names have to be exact. At a minimum the columns the table requires are: bam, TRT, REP, REPrep (often all "1"), Num, and Proportion. Other columns are allows, but are ignored.

bam must match the bam file prefixes (that is the readgroups) of previous steps.  TRT must be C for Controls vs. Z for experiments (if you have other labels use "mutate" and "recode"). Rows with values other than C or Z are ignored. REP is the replicate number and REPrep is a potential technical replicate within that replicate (say a 2nd draw of flies from the same cage).  Num is the number of flies per pool and Proportion the proportion selected. Make sure you change percent to proportions!  Proportion selected should only be associated with "Z" treatments, with "C" -> NA.  Here is an example table, with a few rows to give you an idea.  Note the extra columns associated with this dataset that are ignored, but are part of the table.
```bash
    filesize        bam    A longTRT REP REPrep  Num Proportion TRT
  7953661903 STV1_F_Con STV1     Con   1      1 1205         NA   C
  8479453370 STV1_F_Res STV1     Res   1      1  115     0.0871   Z
  7535619860 STV2_F_Con STV2     Con   2      1 1387         NA   C
  5517579963 STV2_F_Res STV2     Res   2      1  296     0.1540   Z
  8173583358 STV3_F_Con STV3     Con   3      1 1631         NA   C
  6004306536 STV3_F_Res STV3     Res   3      1  174     0.0876   Z
 10790183756 STV4_F_Con STV4     Con   4      1 1628         NA   C
 10309377074 STV4_F_Res STV4     Res   4      1  153     0.0781   Z

```

## Concatenate chromosomes and summarize (10 min)
Up until this point all analyses are done chromosome-by-chromosome for speed.  Now we concatenate and generate some summary figures
```bash
# note the path to the results for each scan above
bash scripts/haps2scan/concat_Chromosome_Scans.Andreas.sh "process/Oct28_24/TEST_A"
bash scripts/haps2scan/concat_Chromosome_Scans.Andreas.sh "process/Oct28_24/TEST_B"
```

## Download the two summary files and some summary plots
```bash
scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/TEST_F/TEST_F.meansBySample.txt .
scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/TEST_F/TEST_F.pseudoscan.txt .
scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/TEST_F/TEST_F.5panel.cM.png .
scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/TEST_F/TEST_F.5panel.Mb.png .
scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/TEST_F/TEST_F.Manhattan.png .

```

The "pseudoscan" and "means" files are pretty rich and likely the only files you need to work with

## Analysis and Visualization Functions

**Note:** The XQTL analysis and plotting functions have been moved to a dedicated R package called [XQTL2.Xplore](https://github.com/tdlong/XQTL2.Xplore).

### Quick Installation

```r
# Install from GitHub with vignettes (important for RStudio users!)
devtools::install_github("tdlong/XQTL2.Xplore", build_vignettes = TRUE)

# Load the package
library(XQTL2.Xplore)
```

### What's Available

The XQTL2.Xplore package provides comprehensive tools for:
- **Genome-wide QTL visualization** with Manhattan plots
- **Regional analysis** with peak refinement tools  
- **Frequency change analysis** across experimental conditions
- **Gene and variant annotation** visualization
- **Publication-ready multi-panel plots**

### Example Usage

```r
# Load example data
data(zinc_hanson_pseudoscan)
data(zinc_hanson_means)

# Genome-wide analysis
XQTL_Manhattan_5panel(zinc_hanson_pseudoscan, cM = FALSE)

# Peak refinement
out <- XQTL_zoom(zinc_hanson_pseudoscan, "chr3R", 18000000, 20000000, 3, 3)

# Publication-ready plot
XQTL_5panel_plot(zinc_hanson_pseudoscan, zinc_hanson_means, 
                 dm6.variants, dm6.ncbiRefSeq.genes, 
                 out$chr, out$start, out$stop)
```

### Documentation

- **Package vignettes** - Complete tutorials and examples
- **Function documentation** - Detailed help for all functions
- **Example datasets** - Ready-to-use data for learning

For the complete analysis workflow and all available functions, please visit the [XQTL2.Xplore repository](https://github.com/tdlong/XQTL2.Xplore).

## Additional Tools

### Sliding Window Heterozygosity Analysis

A useful modification to the pipeline for analyzing low-pass sequence data to assess heterozygosity in strains of interest. This tool is particularly valuable for:

- **Strain validation** - Confirming a tester strain is homozygous before crossing
- **Chromosome extraction verification** - Validating successful chromosome extraction using balancers
- **Quality control** - Assessing strain purity and identifying heterozygous regions

The approach takes BAM files from low-pass samples plus A or B founders (to focus on good SNPs) and generates sliding window plots of heterozygosity. It doesn't require high sequence coverage and has proven useful in practice.

#### Step 1: Create BAM file list

First, create a list of BAM files that includes your tester strains and some inbred founders:

```bash
# Find all BAM files from your samples of interest
find data/bam/June_2025 -name "*.bam" > tester_strains.bams

# Add B-founder BAM files to focus on good SNPs
cat /dfs7/adl/tdlong/fly_pool/XQTL2/helpfiles/founder.bams.txt | grep "B" >> tester_strains.bams

# Run the REFALT generation step
sbatch scripts/raw2bam2REFALT/bam2bcf2REFALT.sh tester_strains.bams process/heterozygosity_analysis
```

#### Step 2: Run heterozygosity analysis

The heterozygosity analysis uses a dedicated script `bam2bcf2REFALT_het.sh` and the `REFALT2HET.awk` utility script (both located in the `scripts/Heterozygosity_tests/` folder).

```bash
# Run the heterozygosity-specific REFALT generation
sbatch scripts/Heterozygosity_tests/bam2bcf2REFALT_het.sh tester_strains.bams process/heterozygosity_analysis
```

This will generate heterozygosity scores for each sample at each SNP position, outputting files like `nHet.chrX.txt`, `nHet.chr2L.txt`, etc.

#### Step 3: Process heterozygosity data and create sliding window summaries

The heterozygosity files are processed using the `process_heterozygosity.R` script to create sliding window summaries across genetic map positions.

```bash
# Process the heterozygosity files
Rscript scripts/Heterozygosity_tests/process_heterozygosity.R process/heterozygosity_analysis/
```

**What this step does:**
- Reads all chromosome-specific heterozygosity files (`nHet.chr*.txt`)
- Adds genetic map positions (cM) using the flymap reference
- Creates sliding window bins (10 bins per chromosome)
- Calculates average heterozygosity per bin per sample
- Outputs a combined table (`het.table.R`) ready for visualization

#### Step 4: Create heterozygosity plots

Generate sliding window heterozygosity plots using the `plot_heterozygosity.R` script. This automatically excludes A and B founder strains and creates publication-ready visualizations.

```bash
# Create heterozygosity plots
Rscript scripts/Heterozygosity_tests/plot_heterozygosity.R process/heterozygosity_analysis/het.table.R
```

**What this step does:**
- Automatically filters out founder strains (A1, A2, B1, B2, etc.)
- Creates chromosome-specific plots with sliding window heterozygosity
- Uses clean, publication-ready styling
- Saves plots as high-resolution PNG files
- Provides sample lists for verification
