# CURRENT STATUS - XQTL2 Exploration Project

## 🚀 **CURRENT STATUS: LIST-FORMAT HAPLOTYPE ESTIMATOR PRODUCTION RUN**

**STATUS**: ✅ **JUICE ANALYSIS COMPLETE** | ✅ **ZINC2 ANALYSIS COMPLETE** | ✅ **SCRIPT REORGANIZATION COMPLETE** | ✅ **PLOTTING INFRASTRUCTURE COMPLETE** | ✅ **SCRIPT PATHS FIXED** | ✅ **SMOOTH_H4 ESTIMATOR INTEGRATED** | 🚀 **LIST-FORMAT HAPLOTYPE ESTIMATOR PRODUCTION RUN**

---

## 📊 **COMPLETED PHASES**

### ✅ **PHASE 1: JUICE DATASET ANALYSIS - COMPLETE**

**Results Summary (9/9 estimators complete):**
```
Method     Parameter  Reliable%  Unreliable%  Failed%
fixed      20kb       72.2%      27.0%        0.9%
fixed      50kb       84.1%      15.6%        0.4%  
fixed      100kb      90.5%      9.5%         0.0%
fixed      200kb      95.6%      4.4%         0.0%
fixed      500kb      100.0%     0.0%         0.0%
adaptive   h4         100.0%     0.0%         0.0%
adaptive   h6         100.0%     0.0%         0.0%
adaptive   h8         100.0%     0.0%         0.0%
adaptive   h10        100.0%     0.0%         0.0%
```

**🎯 KEY ACHIEVEMENT**: Adaptive window algorithm working **perfectly** - 100% reliability across all h_cutoff values!

### ✅ **PHASE 2: JUICE SNP IMPUTATION - COMPLETE**

- ✅ **All 9 parameter combinations completed successfully!**
- ✅ Critical bug fixed and re-run completed
- ✅ All SNP imputation files generated and validated
- ✅ Method evaluation completed with performance metrics

### ✅ **PHASE 3: SCRIPT REORGANIZATION - COMPLETE**

- ✅ **Scripts directory completely reorganized** for clarity and maintainability
- ✅ **Hard-coded JUICE references removed** from all production scripts
- ✅ **Production scripts separated** from debug/testing code
- ✅ **Clean directory structure** with clear documentation

### ✅ **PHASE 4: SCRIPT PATH FIXES - COMPLETE**

- ✅ **Critical path bugs fixed** in all production shell scripts
- ✅ **All scripts now correctly reference** `scripts/production/` directory
- ✅ **Usage examples updated** to show correct paths
- ✅ **Directory structure rules established** and documented

### ✅ **PHASE 5: SMOOTH_H4 ESTIMATOR INTEGRATION - COMPLETE**

- ✅ **New `smooth_h4` estimator created** - 21-position sliding window mean of `adaptive_h4`
- ✅ **Post-processing script developed** (`create_smooth_haplotype_estimator.R`)
- ✅ **SLURM pipeline integration** - Added as 10th array element in parameter table
- ✅ **SNP imputation compatibility** - Updated validation to accept `smooth_h4` format
- ✅ **Quality control implemented** - `estimate_OK` based on 17/21 reliable positions
- ✅ **Normalization applied** - Founder frequencies rescaled to sum to 1.0

---

## 🧬 **SMOOTH_H4 ESTIMATOR - COMPLETE**

### **🎯 Overview**
The `smooth_h4` estimator is a novel post-processing method that applies a 21-position sliding window mean to `adaptive_h4` results, creating smoother haplotype frequency estimates while maintaining biological relevance.

### **🔧 Technical Implementation**
- **Input**: `adaptive_h4` haplotype estimation results
- **Processing**: 21-position sliding window mean applied separately to each founder
- **Quality Control**: `estimate_OK = 1` if ≥17/21 positions are reliable, else `0`
- **Normalization**: Founder frequencies rescaled to sum to 1.0 at each position
- **Output**: `smooth_h4_results_<chr>.RDS` files compatible with SNP imputation

### **📁 Files Created**
```
scripts/production/create_smooth_haplotype_estimator.R  # Post-processing script
helpfiles/production_slurm_params.tsv                   # Updated with 10th element
scripts/production/snp_imputation_from_table.sh         # Updated for smooth_h4
scripts/production/euchromatic_SNP_imputation_single.R  # Updated validation
```

### **🚀 Integration Status**
- ✅ **Parameter table updated**: `chr2R	smooth_h4	4` (10th element)
- ✅ **SLURM array expanded**: `--array=1-10` (was 1-9)
- ✅ **SNP imputation validated**: Accepts `smooth_h4` format
- ✅ **File path handling**: Correctly loads `smooth_h4_results_<chr>.RDS`

### **📊 Expected Benefits**
- **Smoother estimates**: Reduced noise from 21-position averaging
- **Maintained reliability**: Quality control based on underlying `adaptive_h4` data
- **Biological relevance**: Preserves founder frequency relationships
- **Pipeline compatibility**: Seamless integration with existing SNP imputation

### **🎮 Usage Commands**
```bash
# Create smooth_h4 estimator (run once per chromosome)
Rscript scripts/production/create_smooth_haplotype_estimator.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Run SNP imputation for smooth_h4 (10th array element)
sbatch --array=10 scripts/production/snp_imputation_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Run all 10 SNP imputations (including smooth_h4)
sbatch scripts/production/snp_imputation_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
```

---

## ✅ **PHASE 6: LIST-FORMAT HAPLOTYPE ESTIMATOR - COMPLETE**

### **Status**: ✅ **BOTH ADAPTIVE_H4 AND SMOOTH_H4 LIST-FORMAT ESTIMATORS WORKING PERFECTLY**

**What's completed**: 
- **Adaptive_h4 list-format**: ✅ **COMPLETE** - 7,636 rows with correct list-column structure
- **Smooth_h4 list-format**: ✅ **COMPLETE** - 21-position sliding window smoothing working perfectly
- **Output structure**: Perfect list-column format with proper dimensions

**🎯 ACHIEVEMENTS**:
- **Adaptive_h4 list-format**: `process/ZINC2/haplotype_results_list_format/adaptive_window_h4_results_chr2R.RDS` - ✅ **COMPLETE**
- **Smooth_h4 list-format**: `process/ZINC2/haplotype_results_list_format/smooth_h4_results_chr2R.RDS` - ✅ **COMPLETE**
- **List-column structure**: Perfect `<int [8]>`, `<dbl [8]>`, `<dbl [8 × 8]>`, `<chr [8]>` format

**📊 FINAL RESULTS**:
```
# A tibble: 7,636 × 7
   CHROM     pos sample    Groups    Haps      Err           Names    
   <chr>   <dbl> <chr>     <list>    <list>    <list>        <list>   
 1 chr2R 5500000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 2 chr2R 5510000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 3 chr2R 5520000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 4 chr2R 5530000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 5 chr2R 5540000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 6 chr2R 5550000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 7 chr2R 5560000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 8 chr2R 5570000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 9 chr2R 5580000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
10 chr2R 5590000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
# ℹ 7,626 more rows
```

**🎯 KEY TECHNICAL ACHIEVEMENTS**:
1. **Perfect list-column structure**: Each column contains proper vectors/matrices
2. **Correct row count**: 7,636 rows (4 samples × 1,909 positions) after sliding window
3. **Sliding window logic**: Properly handles 21-position windows with quality control
4. **Groups assignment**: `[1,2,3,4,5,6,7,8]` for quality OK, `[1,1,1,1,1,1,1,1]` for quality NOT OK
5. **Haplotype averaging**: Element-wise averaging of frequencies from good positions
6. **Error matrix averaging**: Proper 8×8 matrix averaging with dimension validation
7. **Normalization**: Frequencies properly normalized to sum to 1.0

**📁 FILES CREATED**:
```
scripts/production/create_smooth_haplotype_estimator_list_format.R  # Smooth_h4 list-format script
process/ZINC2/haplotype_results_list_format/
├── adaptive_window_h4_results_chr2R.RDS    # Adaptive_h4 list-format results
└── smooth_h4_results_chr2R.RDS             # Smooth_h4 list-format results
```

**🎮 USAGE COMMANDS**:
```bash
# Create smooth_h4 list-format estimator
Rscript scripts/production/create_smooth_haplotype_estimator_list_format.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Run adaptive_h4 list-format estimator (via SLURM)
sbatch scripts/production/run_list_format_haplotype_estimation_slurm.sh chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R
```

**🎯 NEXT PHASE**: Ready to integrate list-format estimators into complete pipeline

**🔧 TECHNICAL IMPLEMENTATION**:
- **Modified functions**: `scripts/debug/haplotype_estimation_functions_with_groups.R`
- **Key changes**: Added `fulloutput = TRUE` to `lsei` calls to capture error matrices
- **Groups capture**: Direct capture from `cutree` result in clustering step
- **Error matrices**: Capture from `lsei_result$cov` with proper dimension handling
- **Both methods**: Fixed and adaptive methods now return groups and error matrices

**📁 FILES CREATED**:
```
scripts/debug/haplotype_estimation_functions_with_groups.R  # Modified working functions
scripts/debug/run_haplotype_estimation_single_position.R    # Test script for single position
scripts/debug/test_modified_functions.R                     # Test script for new functions
scripts/production/run_haplotype_estimation_list_format.R   # New list-format pipeline
scripts/production/run_snp_imputation_list_format.R         # SNP imputation for list format
scripts/production/run_list_format_pipeline_slurm.sh        # SLURM pipeline for list format
```

**📁 TEST OUTPUT FILES**:
```
process/ZINC2/test_results/
├── adaptive_window_h4_single_position_5400000_results_chr2R.RDS    # Test output for adaptive_h4
├── fixed_window_20kb_single_position_5400000_results_chr2R.RDS     # Test output for fixed_20kb
└── [other test files as needed]
```

**🎯 FILE NAMING CONVENTION**:
- **Production files**: `adaptive_window_h4_results_chr2R.RDS`
- **Test files**: `adaptive_window_h4_single_position_5400000_results_chr2R.RDS`
- **Test directory**: `process/ZINC2/test_results/` (separate from production)
- **Key difference**: Test files include `_single_position_<position>_` to avoid overwriting production

**🧪 TESTING STATUS**:
- **Working function test**: ✅ Successfully captures groups and error matrices
- **List format test**: 🔄 Currently debugging tibble structure
- **Target format**: Working towards correct `[4]` list structure per sample
- **Debug output**: Shows groups being captured correctly from clustering

**🎮 TEST COMMANDS**:
```bash
# Test single position with new functions
Rscript scripts/production/run_haplotype_estimation_list_format.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R 5400000

# Run all positions (production mode)
Rscript scripts/production/run_haplotype_estimation_list_format.R chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R

# Run via SLURM (all chromosomes supported)
sbatch scripts/production/run_list_format_haplotype_estimation_slurm.sh chr2R adaptive 4 process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R

# Test modified functions directly
Rscript scripts/debug/test_modified_functions.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 5400000
```

**📁 OUTPUT MODES**:
- **Single position mode**: Saves to `process/ZINC2/test_results/` with `_single_position_<pos>_` in filename
- **All positions mode**: Saves to `process/ZINC2/haplotype_results_list_format/` with production naming
- **Verbose control**: Add `1` as last argument for debug output
- **Chromosome support**: chr2L, chr2R, chr3L, chr3R, chrX (all 5 chromosomes)

**📁 TEST OUTPUT LOCATIONS**:
```bash
# Test output files are saved to:
process/ZINC2/test_results/adaptive_window_h4_single_position_5400000_results_chr2R.RDS

# Production output files are saved to:
process/ZINC2/haplotype_results/adaptive_window_h4_results_chr2R.RDS

# Key difference: test files include "_single_position_<position>_" in filename
```

**🔍 CURRENT DEBUGGING**:
- **Issue**: Groups being treated as individual founder assignments instead of per-sample vectors
- **Fix**: Modified debug script to group by sample instead of processing each founder separately
- **Target**: 1 row per position with `[4]` lists (one per sample)
- **Progress**: Debug output shows correct groups capture, working on proper tibble structure

**📊 EXPECTED OUTPUT FORMAT**:
```r
# A tibble: 1 × 7
   CHROM    pos      sample       Groups         Haps          Err        Names
   <chr>  <dbl> <list<chr>> <list<list>> <list<list>> <list<list>> <list<list>>
 1 chr2R  5400000        [4]         [4]         [4]         [4]         [4]
```

**🎯 NEXT STEPS**:
1. **Fix tibble structure** to show correct `[4]` format
2. **Validate groups capture** from clustering step
3. **Validate error matrices** from `lsei` with `fulloutput = TRUE`
4. **Test full pipeline** with new list-format estimator
5. **Integrate with SNP imputation** for complete workflow

---

## 🎨 **PLOTTING INFRASTRUCTURE - COMPLETE**

### ✅ **Region Plotting Script** (`plot_summary_region.R`)
- **Three-panel plots**: B1 frequencies, MAE, SNP counts
- **Three methods**: `adaptive_h4`, `fixed_20kb`, `fixed_100kb`
- **Proper NA handling**: Lines break at unreliable estimates (`estimate_OK = FALSE`)
- **Clean styling**: Colorblind-friendly colors, proper y-axis ranges
- **Usage**: `Rscript scripts/production/plot_summary_region.R <chr> <params> <output_dir> <midpoint_10kb>`

### ✅ **Chromosome Plotting Script** (`plot_summary_chromosome.R`)
- **Single method focus**: `adaptive_h4` only
- **Entire chromosome view**: No region filtering
- **Advanced styling**: Transparent black points, LOESS smoothing curves
- **Responsive smoothing**: `span = 0.05` for detailed pattern capture
- **Usage**: `Rscript scripts/production/plot_summary_chromosome.R <chr> <params> <output_dir>`

### 🎯 **Key Plotting Features**:
- **Automatic gap handling**: Lines break at `NA` values (unreliable estimates)
- **Y-axis optimization**: Extended ranges to cover all data points
- **Professional appearance**: Clean themes, appropriate scales
- **Data integrity**: SNP counts from haplotype files, MAE from imputation files

---

## 📋 **FILE FORMAT DOCUMENTATION**

### **🧬 Haplotype Estimation Results - Standard Format**
**File Pattern**: `process/<dataset>/haplotype_results/<method>_results_<chr>.RDS`

**Structure**: Data frame with columns:
- `CHROM`: Chromosome identifier (chr2L, chr2R, chr3L, chr3R, chrX)
- `pos`: Genomic position (integer)
- `sample`: Sample identifier (e.g., Rep01_W_F, Rep02_W_F, etc.)
- `B1_freq`, `B2_freq`, ..., `B8_freq`: Founder haplotype frequencies (0-1, sum to 1.0)
- `estimate_OK`: Boolean flag for reliability (TRUE/FALSE)
- `window_size`: Window size used for estimation (integer, kb)
- `n_snps`: Number of SNPs in window (integer)

**Example**:
```r
# A tibble: 7,636 × 11
   CHROM     pos sample    B1_freq B2_freq B3_freq B4_freq B5_freq B6_freq B7_freq B8_freq estimate_OK window_size n_snps
   <chr>   <dbl> <chr>       <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <lgl>            <dbl>  <int>
 1 chr2R 5500000 Rep01_W_F   0.125   0.125   0.125   0.125   0.125   0.125   0.125   0.125 TRUE                20     45
 2 chr2R 5510000 Rep01_W_F   0.200   0.150   0.100   0.150   0.200   0.100   0.050   0.050 TRUE                20     52
```

### **🧬 Haplotype Estimation Results - List Format**
**File Pattern**: `process/<dataset>/haplotype_results_list_format/<method>_results_<chr>.RDS`

**Structure**: Data frame with list-columns:
- `CHROM`: Chromosome identifier (chr2L, chr2R, chr3L, chr3R, chrX)
- `pos`: Genomic position (integer)
- `sample`: Sample identifier (e.g., Rep01_W_F, Rep02_W_F, etc.)
- `Groups`: List of integers `[8]` - founder group assignments from clustering
- `Haps`: List of doubles `[8]` - founder haplotype frequencies (0-1, sum to 1.0)
- `Err`: List of matrices `[8×8]` - error covariance matrices from lsei
- `Names`: List of strings `[8]` - founder names (B1, B2, ..., B8)

**Example**:
```r
# A tibble: 7,636 × 7
   CHROM     pos sample    Groups    Haps      Err           Names    
   <chr>   <dbl> <chr>     <list>    <list>    <list>        <list>   
 1 chr2R 5500000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
 2 chr2R 5510000 Rep01_W_F <int [8]> <dbl [8]> <dbl [8 × 8]> <chr [8]>
```

**List-Column Access**:
```r
# Access frequencies for first row
frequencies <- data$Haps[[1]]  # Returns vector of 8 frequencies
groups <- data$Groups[[1]]     # Returns vector of 8 group assignments
error_matrix <- data$Err[[1]]  # Returns 8×8 error covariance matrix
founder_names <- data$Names[[1]] # Returns vector of 8 founder names
```

### **🧬 SNP Imputation Results**
**File Pattern**: `process/<dataset>/haplotype_results/snp_imputation_<method>_<chr>.RDS`

**Structure**: Data frame with columns:
- `CHROM`: Chromosome identifier
- `pos`: Genomic position (integer)
- `sample`: Sample identifier
- `snp_id`: SNP identifier
- `ref_allele`: Reference allele
- `alt_allele`: Alternative allele
- `imputed`: Imputed genotype (0, 1, 2, or NA)
- `confidence`: Imputation confidence score (0-1)
- `method`: Imputation method used

### **📊 Summary Files**
**File Pattern**: `process/<dataset>/summary_<chr>_<method>.RDS`

**Structure**: Aggregated results with performance metrics, MAE, coverage statistics, etc.

---

## 📁 **DIRECTORY STRUCTURE & PATH RULES**

### **🌐 Cluster Environment & Relative Paths**
- **ALL DATA AND ANALYSIS HAPPENS ON THE CLUSTER** - never run data analysis commands locally
- All input/output files live on the compute cluster.
- Always use relative paths from the project/process working directory in commands and scripts.
- Do not use absolute local filesystem paths in documentation or usage examples.

### **🚨 CRITICAL WORKFLOW RULE: CLUSTER-ONLY ANALYSIS**
- **NEVER run data analysis commands locally** (Rscript, file operations on data, etc.)
- **ALWAYS follow this workflow**:
  1. Create/edit scripts locally
  2. `git add`, `git commit`, `git push` 
  3. User pulls on cluster and runs analysis
  4. User reports results back
- **Local commands are ONLY for**: git operations, file creation/editing, directory listing
- **Cluster commands are for**: all data analysis, R scripts, file processing
- **This rule prevents**: attempting to access cluster data from local environment

### **🔧 Production Scripts Location**
- **All production scripts**: `scripts/production/`
- **All shell scripts must reference**: `scripts/production/` (not `scripts/`)
- **All R scripts must be called with**: `scripts/production/` prefix

### **📋 Correct Script Paths**:
```bash
# Production shell scripts
scripts/production/haplotype_testing_from_table.sh
scripts/production/snp_imputation_from_table.sh
scripts/production/haplotype_testing_common.sh
scripts/production/haplotype_testing_small_windows.sh

# Production R scripts
scripts/production/run_haplotype_estimation.R
scripts/production/euchromatic_SNP_imputation_single.R
scripts/production/check_snp_imputation_status.R
scripts/production/evaluate_imputation_methods.R
scripts/production/create_summary_file_chunked.R
scripts/production/plot_summary_region.R
scripts/production/plot_summary_chromosome.R
```

### **🚨 Critical Path Rules**:
1. **NEVER use `scripts/`** - always use `scripts/production/`
2. **All shell scripts** must call R scripts with full `scripts/production/` path
3. **Usage examples** in scripts must show correct `scripts/production/` paths
4. **Documentation** must reflect actual script locations

### **📂 Directory Organization**:
```
scripts/
├── production/          # ✅ Production-ready scripts
├── debug/              # 🔍 Debug and testing scripts  
├── archive/            # 📦 Old/backup scripts
├── haps2scan/          # 📊 Haplotype scanning scripts
├── Heterozygosity_tests/ # 🧬 Heterozygosity analysis
└── raw2bam2REFALT/     # 🔄 Data processing scripts
```

### **🎯 Parameter Files**:
- **SLURM parameters**: `helpfiles/production_slurm_params.tsv`
- **JUICE parameters**: `helpfiles/JUICE_haplotype_parameters.R`
- **ZINC2 parameters**: `helpfiles/ZINC2_haplotype_parameters.R`

---

## 🏗️ **SCRIPT REORGANIZATION COMPLETED**

### **New Directory Structure**:
```
scripts/
├── production/          # Main scripts you actually use
│   ├── haplotype_testing_from_table.sh
│   ├── snp_imputation_from_table.sh
│   ├── evaluate_imputation_methods.R
│   ├── check_snp_imputation_status.R
│   ├── create_summary_file_chunked.R
│   ├── plot_summary_region.R
│   └── haplotype_estimation_functions.R
├── debug/              # Debug/testing scripts (one-offs)
├── archive/            # Old/legacy scripts
├── haps2scan/          # Specialized functionality
├── raw2bam2REFALT/     # Data processing
├── Heterozygosity_tests/ # Heterozygosity analysis
└── README.md           # Documentation
```

### **Key Improvements**:
- ✅ **No more hard-coded paths** - everything parameterized
- ✅ **Clear separation** between production and debug code
- ✅ **Easy to find** the scripts you actually need
- ✅ **Documentation** of what each script does
- ✅ **Ready for any dataset** (JUICE, ZINC2, etc.)

---

## 🎯 **ZINC2 ANALYSIS WORKFLOW**

### **Current Pipeline**:
1. ✅ **Haplotype Estimation** - Complete (9 methods)
2. 🔄 **SNP Imputation** - Rerunning with critical bug fix (9 parallel jobs)
3. ⏳ **Method Evaluation** - Waiting for corrected imputation completion
4. ⏳ **Summary Creation** - Waiting for evaluation
5. ⏳ **Visualization** - Waiting for summary

**🔧 Bug Fix Details**:
- **Script**: `scripts/production/euchromatic_SNP_imputation_single.R`
- **Issue**: Not checking `estimate_OK` flag for haplotype reliability
- **Fix**: Now properly flags unreliable estimates as `imputed = NA`
- **Result**: 100% SNP coverage with proper quality control

### **Parameter Combinations** (from `production_slurm_params.tsv`):
```
chr2R	fixed	20
chr2R	fixed	50
chr2R	fixed	100
chr2R	fixed	200
chr2R	fixed	500
chr2R	adaptive	4
chr2R	adaptive	6
chr2R	adaptive	8
chr2R	adaptive	10
chr2R	smooth_h4	4
```

### **Expected Output Structure**:
```
process/ZINC2/
├── RefAlt.chr2R.txt                    # Input SNP data
├── haplotype_results/                   # All results
│   ├── fixed_window_20kb_results_chr2R.RDS
│   ├── fixed_window_50kb_results_chr2R.RDS
│   ├── fixed_window_100kb_results_chr2R.RDS
│   ├── fixed_window_200kb_results_chr2R.RDS
│   ├── fixed_window_500kb_results_chr2R.RDS
│   ├── adaptive_window_h4_results_chr2R.RDS
│   ├── adaptive_window_h6_results_chr2R.RDS
│   ├── adaptive_window_h8_results_chr2R.RDS
│   ├── adaptive_window_h10_results_chr2R.RDS
│   ├── smooth_h4_results_chr2R.RDS
│   ├── snp_imputation_fixed_20kb_chr2R.RDS
│   ├── snp_imputation_fixed_50kb_chr2R.RDS
│   ├── snp_imputation_fixed_100kb_chr2R.RDS
│   ├── snp_imputation_fixed_200kb_chr2R.RDS
│   ├── snp_imputation_fixed_500kb_chr2R.RDS
│   ├── snp_imputation_adaptive_h4_chr2R.RDS
│   ├── snp_imputation_adaptive_h6_chr2R.RDS
│   ├── snp_imputation_adaptive_h8_chr2R.RDS
│   ├── snp_imputation_adaptive_h10_chr2R.RDS
│   └── snp_imputation_smooth_h4_chr2R.RDS
```

---

## 🔑 **KEY ACHIEVEMENTS TO DATE**

1. **✅ Perfect Adaptive Algorithm**: 100% reliability across all parameters
2. **✅ Unified Codebase**: Single architecture for all methods
3. **✅ Robust Testing**: Comprehensive local validation before cluster deployment
4. **✅ Professional Pipeline**: Scalable, well-documented, tidyverse-compliant
5. **✅ Clean Architecture**: Clear separation of concerns, no code duplication
6. **✅ Script Reorganization**: Clean, maintainable script structure
7. **✅ Multi-Dataset Ready**: Pipeline works with any dataset (JUICE, ZINC2, etc.)

---

## ⚠️ **CRITICAL REMINDERS**

### **Testing Protocol ("Laws of Robotics")**
- **Test scripts are source of truth** - production must mimic working tests
- **Always test locally first** - cluster debugging is expensive
- **Never overwrite production data** - use separate test outputs

### **Adaptive Algorithm Sensitivity**
- Complex algorithm with constraint accumulation
- Test thoroughly before any changes
- Verify different h_cutoff values produce different results

### **Script Organization**
- **Use scripts from `production/` directory** for all analysis
- **Debug scripts are in `debug/`** - don't use for production
- **All paths are now parameterized** - no hard-coded references

---

## 📝 **ZINC2 ANALYSIS COMMANDS**

### **Current Pipeline**:
```bash
# Haplotype testing (COMPLETED)
sbatch scripts/production/haplotype_testing_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# SNP imputation (CURRENTLY RUNNING)
sbatch scripts/production/snp_imputation_from_table.sh helpfiles/production_slurm_params.tsv helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
```

### **Next Steps** (after SNP imputation completes):
```bash
# Check completion status
Rscript production/check_snp_imputation_status.R helpfiles/production_slurm_params.tsv process/ZINC2

# Evaluate methods
Rscript production/evaluate_imputation_methods.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Create summary
Rscript production/create_summary_file_chunked.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Plot specific region
Rscript production/plot_summary_region.R chr2R helpfiles/ZINC2_haplotype_parameters.R process/ZINC2 870
```

### **Monitoring Commands**:
```bash
# Check job status
squeue -u $USER

# Check output directory
ls -la process/ZINC2/haplotype_results/

# Check specific file completion
Rscript scripts/production/check_snp_imputation_status.R helpfiles/production_slurm_params.tsv process/ZINC2
```

---

## 🧹 **PROJECT ORGANIZATION STATUS**

### **Scripts Directory**: ✅ **COMPLETELY REORGANIZED**
- **Production scripts**: Clean, parameterized, ready for any dataset
- **Debug scripts**: Separated but accessible for troubleshooting
- **Archive scripts**: Preserved for reference
- **Documentation**: Clear README with usage examples

### **Data Organization**: ✅ **CLEAN AND SCALABLE**
- **JUICE dataset**: Complete analysis in `process/JUICE/`
- **ZINC2 dataset**: Active analysis in `process/ZINC2/`
- **Parameter files**: Dataset-specific configurations
- **Output structure**: Consistent across all datasets

---

*Last Updated: 2025-01-19 - List-Format Haplotype Estimators Complete, Perfect List-Column Structure Achieved, Smooth_H4 Sliding Window Working, Ready for Pipeline Integration*
