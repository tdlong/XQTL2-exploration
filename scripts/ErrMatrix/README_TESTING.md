# Wide Format Optimization Testing

## 🎉 FINAL STATUS: 19.5x Speedup Achieved!

**See `FINAL_STATUS.md` for comprehensive documentation and ready-to-use commands.**

## Current Status (Updated)

**BASE_VAR_WIDE ESTIMATOR TIMING COMPLETED** ✅🚀
- **Script**: `BASE_VAR_WIDE.R` 
- **Command**: `Rscript scripts/ErrMatrix/BASE_VAR_WIDE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose`
- **Results**: 53.19 seconds (0.89 minutes) for 3000 function calls = 0.0177 seconds per call
- **Status**: Successfully completed on cluster
- **Performance**: **19.5x FASTER** than BASE.R (0.3448 sec/call vs 0.0177 sec/call)!

**BASE_VAR ESTIMATOR TIMING COMPLETED** ✅
- **Script**: `BASE_VAR.R` 
- **Command**: `Rscript scripts/ErrMatrix/BASE_VAR.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose`
- **Results**: 1597.4 seconds (26.62 minutes) for 3000 function calls = 0.5325 seconds per call
- **Status**: Successfully completed on cluster
- **Comparison**: 54% slower than BASE.R (0.3448 sec/call vs 0.5325 sec/call)

**BASE ESTIMATOR TIMING COMPLETED** ✅
- Using `BASE.R` with command: `Rscript scripts/ErrMatrix/BASE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose`
- Testing the BASE estimator performance with full adaptive window algorithm
- Using existing working code (no rewrites) - principle: don't rewrite working code for benchmarking
- **RESULTS**: 1034.53 seconds (17.24 minutes) for 3000 function calls = 0.3448 seconds per call

**LOCAL TESTING COMPLETED**: Haplotype error workbench validation
- Successfully tested `run_100_with_dataframe(h_cutoff = 4)` locally
- **100% convergence rate**: All 100 simulations reached 8 groups (all founders distinguishable)
- **Proper variance/covariance estimation**: All error matrices well-conditioned (finite kappa values)
- **Performance metrics**: Mean haplotype error = 0.85%, mean log10(kappa) = 2.02
- **Algorithm efficiency**: Mean 2.0 window transitions to reach full separation
- **Confirmed working**: `estimate_haplotypes_list_format_sim` is the current, validated haplotype estimator

**NEW FUNCTION ADDED**: `est_haps_var` - Advanced haplotype estimation with variance/covariance
- **Purpose**: Enhanced haplotype estimation with sophisticated error matrix estimation
- **Key Features**: 
  - Uses genomic distance-based windowing (like EHLF) instead of SNP count-based
  - Progressive V matrix construction with pending covariance handling
  - Pooled covariance estimation for grouped founders
  - Constraint accumulation across window sizes
- **Window Sizes**: 10kb, 20kb, 50kb, 100kb, 200kb, 500kb (genomic distances)
- **Windowing Logic**: testing_position ± window_size/2 (centered on test position)
- **Status**: Ready for testing with real genomic data

### Reshaping/Output Integrity (Updated)
- The prior mismatch between unreshaped adaptive `Err` matrices and the reshaped `adapt_h4` output has been fixed inside `BASE_VAR_WIDE.R`’s smoothing/reshaping flow.
- The fix preserves per-sample alignment before reshaping and copies the full 8×8 `Err` matrix without recomputation.
- First-row strict check shows perfect agreement (names, haps, full `Err` matrix, and trace).
- Comparison scripts hardened to never abort on problematic rows; they emit NA metrics and CSVs for inspection.

## Cluster Deployment Process

**CRITICAL CONSTRAINT**: All testing and execution happens on the cluster, not locally!

**DO NOT** attempt to run R scripts locally using `run_terminal_cmd` - they must be executed on the high-performance cluster.

**Local environment is only for**:
- Code editing and file management
- Git operations (add, commit, push)
- Documentation updates

### To Deploy Changes:
1. **Commit and push changes locally**:
   ```bash
   git add .
   git commit -m "Description of changes"
   git push origin main
   ```

2. **Run on cluster**:
   - SSH into cluster
   - Pull latest changes: `git pull origin main`
   - Execute the desired testing script

### Current Testing Location
- **Cluster path**: `/path/to/cluster/XQTL2_exploration/scripts/ErrMatrix/`
- **Active script**: `BASE.R` (existing working code)
- **Status**: Currently running BASE estimator timing tests

## Where to Run From
Run all commands from the **project root directory** (`/Users/tdlong/Desktop/Cursor/XQTL2_exploration/`):

```bash
cd /Users/tdlong/Desktop/Cursor/XQTL2_exploration
```

## Testing Scripts

### 1. Test Production Version
```bash
Rscript scripts/ErrMatrix/test_production.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

### 2. Test Wide Format Version  
```bash
Rscript scripts/ErrMatrix/test_wide_format.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

### 3. Benchmark Both Versions
```bash
Rscript scripts/ErrMatrix/benchmark_wide_vs_production.R chr2R process/JUICE helpfiles/JUICE_haplotype_parameters.R
```

### 4. BASE Estimator Timing (Currently Running)
```bash
Rscript scripts/ErrMatrix/BASE.R chr2R adaptive 4 process/JUICE helpfiles/JUICE_haplotype_parameters.R --debug --nonverbose
```

### 5. Haplotype Error Workbench Testing (Local - COMPLETED)
```bash
Rscript -e "source('scripts/ErrMatrix/haplotype_error_workbench.R'); run_100_with_dataframe(h_cutoff = 4)"
```

### 6. est_haps_var Function Testing (Ready for Testing)
```bash
# Test with simulated data (pos = -99)
Rscript -e "source('scripts/ErrMatrix/haplotype_error_workbench.R'); run_100_with_dataframe(h_cutoff = 4)"

# Test with real genomic data (testing_position = actual genomic position)
# Note: Requires real df3 data with genomic coordinates
```

## What Each Script Does

### Current Production Script
- **`BASE_VAR_WIDE.R`**: Current production script with 19.5x speedup (0.0177 sec/call)
  - **WIDE FORMAT OPTIMIZATION**: Eliminates long↔wide pivoting overhead
  - **Data Structure**: POS, founder1, founder2, ..., foundern, sample1, sample2, ..., sampleM
  - **Performance**: Direct column access instead of pivot_longer/pivot_wider
  - **Complete Workflow**: Adaptive estimation + smoothing + output formatting
  - **Reshaping fix included here**: The corrected reshaping logic is implemented inside this script. Re-running this single script (via SLURM wrapper) produces correct `adapt_h4` outputs.
  - **SLURM Ready**: Use `run_all_chroms.slurm` for cluster deployment

### Archived Scripts (in `archive/` folder)
- **`BASE.R`**: Original working BASE estimator (timing completed: 0.3448 sec/call)
- **`BASE_VAR.R`**: Advanced variance/covariance estimator (timing completed: 0.5325 sec/call)

### Supporting Scripts
- **`haplotype_error_workbench.R`**: Main simulation and testing framework
  - **`est_haps_var`**: Advanced haplotype estimation with genomic distance-based windowing and progressive variance/covariance estimation
  - **`run_100_with_dataframe`**: Wrapper for 100 simulations with proper validation
- **`run_all_chroms.slurm`**: SLURM wrapper for cluster deployment across all chromosomes

## Safe Testing
- All scripts use debug mode (limited data)
- Wide format saves to `haplotype_results_list_format_wide_optimized/` (won't overwrite production)
- No risk to existing production data
- BASE.R uses existing working code (no rewrites) - principle: don't rewrite working code for benchmarking

## Expected Output
- Timing information for each version
- Speedup measurement (hopefully 5-10x faster)
- Verification that results are similar between versions
- **BASE estimator performance metrics (COMPLETED)**: 0.3448 seconds per call for 3000 function calls
- **BASE_VAR estimator performance metrics (COMPLETED)**: 0.5325 seconds per call for 3000 function calls (54% slower than BASE)
- **BASE_VAR_WIDE estimator performance metrics (COMPLETED)**: 0.0177 seconds per call for 3000 function calls (**19.5x FASTER than BASE!**)
- Haplotype estimator validation results (100% convergence, proper var/cov estimation)

## Recent Accomplishments

### 1. Function Development & Integration
- ✅ **Created `est_haps_var`**: Advanced haplotype estimator with genomic distance-based windowing
- ✅ **Integrated into BASE_VAR.R**: Replaced original `estimate_haplotypes_list_format` with `est_haps_var`
- ✅ **Maintained compatibility**: Same input/output interface as original function
- ✅ **Removed delegation logic**: Always runs advanced estimation (no fallback to production)

### 2. Code Architecture Improvements
- ✅ **Genomic distance-based windowing**: Uses `testing_position ± window_size/2` like EHLF
- ✅ **EHLF window sizes**: 10kb, 20kb, 50kb, 100kb, 200kb, 500kb
- ✅ **Advanced variance/covariance**: Progressive V matrix construction and pooled covariance
- ✅ **Constraint accumulation**: Builds constraints across window sizes for better estimation

### 3. Speed Optimization (COMPLETED) 🚀
- ✅ **Wide Format Optimization**: Created `BASE_VAR_WIDE.R` to eliminate pivoting overhead
- ✅ **Data Structure**: POS, founder1, founder2, ..., foundern, sample1, sample2, ..., sampleM
- ✅ **Direct Access**: No more `pivot_longer`/`pivot_wider` - direct column access
- ✅ **MASSIVE SPEEDUP**: **19.5x FASTER** than BASE.R (0.3448 sec/call → 0.0177 sec/call)
- ✅ **Pre-subset Strategy**: Filter data once per position instead of 6 times per position×sample

### 4. Testing & Validation
- ✅ **Local validation**: `run_100_with_dataframe` with 100% convergence and proper var/cov estimation
- ✅ **BASE.R benchmarking**: Completed (0.3448 seconds per call for 3000 function calls)
- ✅ **BASE_VAR.R benchmarking**: Completed (0.5325 seconds per call for 3000 function calls)
- ✅ **BASE_VAR_WIDE.R benchmarking**: Completed (0.0177 seconds per call for 3000 function calls) - **19.5x FASTER!**

## Next Steps

### Immediate (Performance Analysis) ✅
1. **✅ Compare performance**: BASE.R (0.3448 sec/call) vs BASE_VAR.R (0.5325 sec/call) - 54% slower
2. **✅ Analyze trade-offs**: Advanced variance/covariance estimation vs performance cost
3. **✅ Test BASE_VAR_WIDE.R**: Wide format optimization delivered **19.5x speedup!**
4. **✅ Document results**: Updated README with performance comparison

### Performance Summary 🚀
- **BASE.R**: 0.3448 sec/call (baseline)
- **BASE_VAR.R**: 0.5325 sec/call (54% slower - advanced var/cov estimation)
- **BASE_VAR_WIDE.R**: 0.0177 sec/call (**19.5x FASTER** - wide format + pre-subset optimization)

### Production Deployment
1. **SLURM Wrapper**: Created `run_all_chroms.slurm` for cluster deployment
   - Processes all 5 chromosomes (chrX, chr2L, chr2R, chr3L, chr3R)
   - Configurable via command line arguments: `<param_file> <output_dir> [parameter]`
   - Runs BASE_VAR_WIDE.R in production mode (non-debug, non-verbose)
   - Complete workflow: adaptive estimation + smoothing + output formatting

2. **Command**: 
   ```bash
   # JUICE parameters (default parameter=4)
   sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/JUICE_haplotype_parameters.R process/JUICE
   
   # ZINC2 parameters (default parameter=4)
   sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/ZINC2_haplotype_parameters.R process/ZINC2
   
   # Custom parameter value
   sbatch scripts/ErrMatrix/run_all_chroms.slurm helpfiles/JUICE_haplotype_parameters.R process/JUICE 6
   ```
   - Creates production-ready output files in specified directory
   - Generates both adaptive and smooth results in proper format
   - Expected runtime: ~1-2 hours per chromosome (19.5x faster than BASE.R)

### Sanity-Check & Validation Tools (New)
- Summary metrics across all rows:
  ```bash
  Rscript scripts/ErrMatrix/compare_adaptive_vs_reshaped.R <chr> <output_dir>
  ```
- Strict equality + CSV of mismatches:
  ```bash
  Rscript scripts/ErrMatrix/strict_compare_adapt_vs_reshaped.R <chr> <output_dir> 1e-10
  ```
Both scripts are hardened and will not abort on malformed rows.

### Future Development
1. **Integration**: Consider integrating best features into production pipeline
2. **Validation**: Test with different datasets and parameters
3. **Monitoring**: Track performance across all chromosomes
