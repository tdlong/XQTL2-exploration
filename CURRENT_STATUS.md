# CURRENT STATUS - XQTL2 Exploration Project

## Current Phase: Phase 4 (Method evaluation and comparison) + Debugging Adaptive Window Algorithm

## Recent Progress (Last 24 hours)

### ✅ Completed Tasks
1. **Full pipeline execution** - Haplotype estimation and SNP imputation completed for chr2R
2. **Comprehensive method evaluation** - All 12 estimators (6 fixed + 6 adaptive) evaluated
3. **Sliding window analysis integration** - Regional performance analysis (250 SNP windows, 50 SNP steps)
4. **Performance optimization** - Sliding window calculation optimized with vectorized tidyverse operations
5. **Documentation updates** - ADAPTIVE_WINDOW_ALGORITHM.md created and updated

### 🔍 Current Focus: Debugging Adaptive Window Algorithm

**The Problem**: All adaptive window methods (h4, h6, h8, h10) are producing identical haplotype estimates and thus identical SNP imputation performance.

**Initial Investigation**: 
- Created `scripts/check_adaptive_estimates.R` to compare results
- Created `scripts/analyze_adaptive_results.R` to analyze cluster results
- Cluster analysis confirmed: identical founder frequencies, always 8 groups, only 2429 positions processed

**Code Analysis**:
- Examined git history and current production code
- **DISCOVERY**: Current production version DOES have hierarchical clustering algorithm
- **MYSTERY**: Why are cluster results identical despite correct algorithm?

**Recent Breakthrough**:
- Created `scripts/test_working_adaptive.R` with verbose output
- **CRITICAL BUG IDENTIFIED**: Redundant constraint accumulation
- **The Problem**: 10kb window (1 founder group) adds constraint "sum=1" which is already enforced
- **The Impact**: Matrix singularity, LSEI failures, unexpected behavior

**Algorithm Understanding Clarified**:
- Start small (10kb) and grow progressively
- Only run LSEI when groups meaningfully change (not just count, but composition)
- Skip LSEI when groups unchanged, reuse previous constraints
- Handle tree reshuffling gracefully (tree structure may change with more SNPs)
- Continue until full founder separation OR max window size (500kb)

## Key Findings

### Adaptive Window Algorithm
- **IS working correctly** in current production code
- **Progressive window expansion** with hierarchical clustering
- **Constraint accumulation** from smaller windows
- **Bug identified**: Redundant constraint accumulation causing matrix issues

### Performance Results
- **Fixed windows**: 10kb performs best, larger windows degrade performance
- **Adaptive windows**: Superior to fixed windows across all metrics
- **All adaptive methods identical**: Suggests implementation bug, not algorithm flaw

## Next Steps

### Immediate (Next 2-4 hours)
1. **Fix redundant constraint bug** in test script
2. **Test fixed algorithm** on same positions with h4 vs h10
3. **Verify different results** for different h_cutoff values
4. **Update production code** if bug confirmed

### Short Term (Next 1-2 days)
1. **Re-run adaptive window jobs** on cluster with fixed code
2. **Verify different performance** across h_cutoff values
3. **Extend to other chromosomes** (chr2L, chr3L, chr3R, chrX)
4. **Genome-wide analysis** of all methods

### Medium Term (Next week)
1. **Parameter optimization** - find optimal h_cutoff values
2. **Method comparison** - comprehensive evaluation across genome
3. **Performance profiling** - runtime and accuracy trade-offs
4. **Documentation** - final algorithm description and usage guide

## Technical Details

### Files Modified Today
- `scripts/test_working_adaptive.R` - Verbose test script with redundant constraint bug
- `ADAPTIVE_WINDOW_ALGORITHM.md` - Comprehensive algorithm documentation
- `CURRENT_STATUS.md` - This file, updated with current debugging status

### Key Commands & Pipeline Execution

#### Current Debugging
```bash
# Test the working adaptive algorithm with verbose output
git pull origin main
Rscript scripts/test_working_adaptive.R
```

#### Full Pipeline (when ready)
```bash
# 1. Haplotype estimation
sbatch scripts/haplotype_testing_from_table.sh

# 2. SNP imputation (optional)
# Edit haplotype_testing_from_table.sh to include SNP imputation

# 3. Evaluation
Rscript scripts/evaluate_haplotype_methods.R
```

### Parameter Table Format
```tsv
chr2R	fixed	10
chr2R	fixed	25
chr2R	fixed	50
chr2R	adaptive	4
chr2R	adaptive	6
chr2R	adaptive	8
chr2R	adaptive	10
```

## Workflow Constraints

### Cluster Access
- **Results stay on cluster** - cannot download large .RDS files locally
- **Local development** - scripts developed locally, tested on cluster
- **Version control** - git push/pull for code synchronization

### Data Processing
- **Raw REFALT files** - direct processing, no intermediate files
- **Euchromatin boundaries** - chr2R: 5,398,184 to 24,684,540 bp
- **10kb step intervals** - comprehensive chromosome coverage

## Key Commands & Pipeline Execution

### 1. Submit Haplotype Estimation Jobs
```bash
# Submit array job for all methods
sbatch scripts/haplotype_testing_from_table.sh

# Monitor progress
squeue -u $USER
squeue -u $USER --array

# Check specific job
squeue -j <job_id>
```

### 2. Check Job Output
```bash
# Check output files
ls -la process/JUICE/adaptive_window_h*_results_chr2R.RDS
ls -la process/JUICE/fixed_window_*_results_chr2R.RDS

# Monitor specific job output
tail -f process/JUICE/haplotype_pipeline_<job_id>.out
tail -f process/JUICE/haplotype_pipeline_<job_id>.err
```

### 3. Run SNP Imputation (Optional)
```bash
# Edit haplotype_testing_from_table.sh to include SNP imputation
# Or run manually for specific methods
Rscript scripts/euchromatic_SNP_imputation_single.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE adaptive 4
```

### 4. Evaluate Results
```bash
# Comprehensive evaluation
Rscript scripts/evaluate_haplotype_methods.R

# Check output
ls -la process/JUICE/evaluation_*
```

### 5. Individual Script Execution
```bash
# Test specific method
Rscript scripts/REFALT2haps.FixedWindow.Single.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE 10

Rscript scripts/REFALT2haps.AdaptWindow.Single.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE 4
```

### 6. Monitoring and Debugging
```bash
# Check job status
sacct -j <job_id> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,MaxVMSize

# Cancel specific jobs
scancel <job_id>

# Cancel all array jobs
scancel -n haplotype_pipeline
```

### 7. Small-Scale Testing
```bash
# Test on subset of data
Rscript scripts/test_working_adaptive.R

# Debug specific issues
Rscript scripts/check_adaptive_estimates.R
```

### 8. File Structure
```
process/JUICE/
├── RefAlt.chr2R.txt          # Raw REFALT data
├── adaptive_window_h*_results_chr2R.RDS  # Adaptive window results
├── fixed_window_*_results_chr2R.RDS      # Fixed window results
├── snp_imputation_*_chr2R.RDS            # SNP imputation results
└── evaluation_*_chr2R.RDS                # Evaluation results
```

## Recent Debugging Insights

### The Redundant Constraint Bug
**What Happens**:
1. 10kb window: All founders in 1 group → adds constraint "sum=1"
2. This constraint is **redundant** with the base constraint already enforced
3. **Matrix singularity** occurs due to colinear rows
4. LSEI may fail or behave unexpectedly

**Why This Matters**:
- Could explain why cluster results are identical
- Matrix issues might cause algorithm to fall back to default behavior
- Different h_cutoff values might all hit the same bug

**The Fix**:
- Only accumulate constraints when groups meaningfully change
- Skip LSEI when groups unchanged
- Detect meaningful changes by group composition, not just count

## Current Questions

1. **Is the redundant constraint bug the root cause** of identical cluster results?
2. **Will fixing this bug** restore different performance across h_cutoff values?
3. **Should we update production code** or just test the fix locally first?
4. **What's the best approach** for detecting meaningful group changes?

## Success Metrics

### Algorithm Working Correctly
- [ ] Different h_cutoff values produce different results
- [ ] Progressive window expansion shows meaningful group changes
- [ ] Constraint accumulation works without matrix issues
- [ ] Full founder separation achieved or max window size reached

### Pipeline Performance
- [ ] All methods complete successfully
- [ ] Different performance metrics across h_cutoff values
- [ ] Reasonable runtime (not stuck in "doom loops")
- [ ] Full chromosome coverage achieved

### Code Quality
- [ ] No redundant constraints
- [ ] Smart LSEI execution (only when needed)
- [ ] Proper error handling and edge cases
- [ ] Clear documentation and testing
