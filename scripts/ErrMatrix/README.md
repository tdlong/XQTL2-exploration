# ErrMatrix Directory - Production Haplotype Estimation

This directory contains the production-ready haplotype estimation pipeline with optimized performance and bug fixes.

## Production Files

### BASE_VAR_WIDE.R
**Purpose**: Production-ready haplotype estimation pipeline with optimized performance  
**Performance**: 19.5x faster than original BASE.R (0.0177 seconds per call)  
**Description**: Complete workflow with adaptive estimation, smoothing, and output formatting  
**Usage**: 
```bash
# Full pipeline for all chromosomes
sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/ZINC2_haplotype_parameters.R process/ZINC2

# Single chromosome testing
Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose
```

**Key Features**:
- ✅ **Optimized Performance**: 19.5x speedup over baseline
- ✅ **Bug-Free Smoothing**: Correctly produces exactly 20 fewer positions
- ✅ **Data Integrity**: Reshaping preserves exact values without recomputation
- ✅ **Production Ready**: Complete workflow with all output formats

**Output Files**:
- `adaptive_window_h4_results_chr*.RDS`: Original adaptive results
- `smooth_h4_results_chr*.RDS`: Smoothed results (20 fewer positions)
- `smooth_h4/R.haps.chr*.out.rds`: Reshaped smoothed data
- `adapt_h4/R.haps.chr*.out.rds`: Reshaped adaptive data

### Supporting Scripts

#### run_smoothing_only.R
**Purpose**: Test smoothing fixes without full pipeline rerun  
**Usage**: 
```bash
Rscript scripts/ErrMatrix/run_smoothing_only.R chr3R process/ZINC2 helpfiles/ZINC2_haplotype_parameters.R
```

#### compare_adaptive_vs_reshaped.R
**Purpose**: Verify data integrity between original and reshaped formats  
**Usage**: 
```bash
Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R chr3R process/ZINC2
```

#### run_all_chroms.slurm
**Purpose**: SLURM job script for running all chromosomes on cluster  
**Usage**: 
```bash
sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
```

## Bug Fixes Applied
- ✅ **Smoothing Position Count**: Fixed incorrect 40-position reduction (now correctly 20 positions)
- ✅ **Window Logic**: Corrected to only process positions 11 through (n-10) for 21-position window
- ✅ **Data Integrity**: Verified reshaping preserves exact values without recomputation

## Performance Achievements
- **19.5x speedup** over original implementation
- **0.0177 seconds per call** (vs 0.3448 for original)
- **Complete workflow** in single optimized script
- **Production-ready** with all required output formats
