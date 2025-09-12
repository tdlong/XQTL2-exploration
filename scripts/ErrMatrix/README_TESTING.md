# Wide Format Optimization Testing

## Current Status (Updated)

**BASE ESTIMATOR TIMING COMPLETED**: Results available
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

## Cluster Deployment Process

**IMPORTANT**: All testing is done on the cluster, not locally.

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

- **test_production.R**: Runs production version in debug mode (100 positions × 1 sample)
- **test_wide_format.R**: Runs wide format version in debug mode (100 positions × 1 sample)  
- **benchmark_wide_vs_production.R**: Runs both versions on same data and compares performance
- **BASE.R**: Existing working code for BASE estimator timing (completed: 0.3448 sec/call)
- **BASE_VAR.R**: Clean duplicate of BASE.R for variance/covariance estimation work
- **haplotype_error_workbench.R**: Validated haplotype estimator with proper var/cov estimation (tested locally)
  - **`est_haps_var`**: Advanced haplotype estimation with genomic distance-based windowing and progressive variance/covariance estimation

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
- Haplotype estimator validation results (100% convergence, proper var/cov estimation)
