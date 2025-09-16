# Quick Start Guide - BASE_VAR_WIDE.R

## ⚠️ CRITICAL INFORMATION - READ FIRST
- **FILES ARE BIG**: Output files can be 50-100MB+ per chromosome
- **RUNS TAKE A LONG TIME**: Full pipeline takes several hours to days
- **ALL DATA IS ON CLUSTER**: You can only commit and push changes locally
- **ASK BEFORE OVERWRITING**: Files take days to create - always ask before overwriting!
- **CLUSTER ACCOUNT REQUIRED**: You MUST specify account and partition, otherwise jobs default to limited "free" resources!

## 🚀 Ready to Use Commands

### Cluster Deployment (All Chromosomes)
```bash
# JUICE dataset (default)
sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/JUICE_haplotype_parameters.R process/JUICE

# ZINC2 dataset
sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Custom parameter
sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/JUICE_haplotype_parameters.R process/JUICE 6
```

### Local Testing (Single Chromosome)
```bash
# Debug mode (500 positions × 4 samples)
Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose

# Production mode (all positions)
Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --nonverbose
```

### 🔧 Testing Fixes (Smoothing Only)
**Use this to test smoothing fixes without rerunning the entire pipeline!**
```bash
# Test smoothing fix on existing adaptive results
Rscript scripts/ErrMatrix/run_smoothing_only.R chr3R process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R

# This will:
# 1. Load existing adaptive_window_h4_results_chr3R.RDS
# 2. Run only the smoothing step with the fix
# 3. Generate the 3 output files
# 4. Verify position counts are correct (should be 20 fewer positions)
```

## 📊 Performance
- **19.5x faster** than original BASE.R
- **0.0177 seconds per call** (vs 0.3448 for BASE.R)
- **Complete workflow**: Adaptive estimation + smoothing + output formatting

## 📁 Data Locations
**ALL DATA IS ON CLUSTER - You can only commit and push changes locally!**

### ZINC2 Dataset
```
process/ZINC2/haplotype_results_list_format/
├── adaptive_window_h4_results_chr3R.RDS  (original adaptive results)
├── smooth_h4_results_chr3R.RDS           (smoothed results)
├── smooth_h4/R.haps.chr3R.out.rds       (reshaped smoothed)
└── adapt_h4/R.haps.chr3R.out.rds        (reshaped adaptive)
```

### JUICE Dataset
```
process/JUICE/haplotype_results_list_format/
├── adaptive_window_h4_results_chr2R.RDS
├── smooth_h4_results_chr2R.RDS
├── smooth_h4/R.haps.chr2R.out.rds
└── adapt_h4/R.haps.chr2R.out.rds
```

### ⚠️ File Safety
- **Files are 50-100MB+ each** - takes days to regenerate
- **Always ask before overwriting** existing results
- **Use `--only-smoothing` mode** to test fixes without full rerun
- **Backup important files** before making changes

## 📋 What It Does
1. **Adaptive Estimation**: Uses 6 window sizes (10kb-500kb) with genomic distance filtering
2. **Variance/Covariance**: Advanced error matrix estimation with progressive V matrix
3. **Smoothing**: 21-position sliding window with quality-based averaging
4. **Output**: Creates production-ready files in expected format

## 🔍 Monitor Progress
```bash
# Check job status
squeue -u $USER

# View logs
tail -f logs/base_var_wide_all_*.out
```

## ✅ Verify Results (After Completion)
```bash
# Run comparison for all chromosomes
sbatch scripts/ErrMatrix/run_comparison_all_chroms.slurm process/ZINC2

# Or test with a limit first
sbatch scripts/ErrMatrix/run_comparison_all_chroms.slurm process/ZINC2 100

# Check comparison results
tail -f logs/compare_adapt_reshaped_*.out
```

**Ready to go!** 🎉
