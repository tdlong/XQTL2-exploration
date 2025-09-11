# ErrMatrix Directory - Haplotype Estimation Optimization

This directory contains experimental code for optimizing haplotype estimation performance.

## Files

### BASE.R
**Purpose**: Baseline production code for timing comparison  
**Performance**: 0.5562 seconds per function call (100 calls = 55.62 seconds)  
**Description**: Exact copy of production `complete_haplotype_workflow.R` with smoothing skipped  
**Usage**: 
```bash
Rscript scripts/ErrMatrix/BASE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose
```

**Key Functions**:
- `estimate_haplotypes_list_format()`: Core adaptive window haplotype estimation
- `process_refalt_data()`: Load and process RefAlt.txt files
- `run_adaptive_estimation()`: Orchestrate adaptive estimation workflow

**Output**: 
- `adaptive_window_h4_results_chr2R.RDS`: Results in production format
- Timing information printed to console

### BASE_VAR.R
**Purpose**: Production code with progressive error matrix function swapped in  
**Performance**: 0.5398 seconds per function call (100 calls = 53.98 seconds)  
**Description**: Same as BASE but uses the progressive error matrix algorithm for correct var/covariance estimation  
**Usage**: 
```bash
Rscript scripts/ErrMatrix/BASE_VAR.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose
```

**Key Functions**:
- `estimate_haplotypes_list_format()`: Progressive error matrix version (delegates to production for real data)
- `estimate_haplotypes_list_format_prod()`: Production function for delegation
- `process_refalt_data()`: Load and process RefAlt.txt files
- `run_adaptive_estimation()`: Orchestrate adaptive estimation workflow

**Output**: 
- `adaptive_window_h4_results_chr2R.RDS`: Results in production format
- Timing information printed to console
- **Improved error matrix estimation** with progressive var/covariance handling

### Optimization Goal
Target: 5-10x speedup by eliminating expensive `pivot_wider` and `arrange` operations in the core haplotype estimation function.

**Current Status**: 
- ✅ Baseline established: 0.5562 seconds per call (BASE.R)
- ✅ Progressive error matrix implemented: 0.5398 seconds per call (BASE_VAR.R)
- 🔄 Next: Create optimized version with pre-conditioned wide-format data
- 🔄 Next: Benchmark optimized vs baseline

### Performance Bottlenecks Identified
1. **`pivot_wider`** in `estimate_haplotypes_list_format()` - called for every position×sample
2. **`arrange`** operations in data processing
3. **Repeated data reshaping** for each window size

### Optimization Strategy
1. Pre-process entire chromosome into wide format once
2. Subset windows from pre-conditioned data
3. Eliminate redundant `pivot_wider`/`arrange` operations
4. Maintain exact same input/output contract as production code
