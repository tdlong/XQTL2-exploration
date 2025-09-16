# Quick Start Guide - BASE_VAR_WIDE.R

## 🚀 Ready to Use Commands

### ⚠️ IMPORTANT: Cluster Account Settings
**On this cluster, you MUST specify account and partition, otherwise jobs default to limited "free" resources!**

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

## 📊 Performance
- **19.5x faster** than original BASE.R
- **0.0177 seconds per call** (vs 0.3448 for BASE.R)
- **Complete workflow**: Adaptive estimation + smoothing + output formatting

## 📁 Output Files
```
process/JUICE/haplotype_results_list_format/
├── adaptive_window_h4_results_chr2R.RDS
├── smooth_h4_results_chr2R.RDS
├── smooth_h4/R.haps.chr2R.out.rds
└── adapt_h4/R.haps.chr2R.out.rds
```

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
