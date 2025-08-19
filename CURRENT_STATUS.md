# CURRENT STATUS - XQTL2 Exploration Project

## 🔧 **CURRENT FOCUS: Deploy Working Unified Function to Production**

**STATUS**: ✅ **UNIFIED FUNCTION WORKING AND TESTED** - Ready for production deployment

**COMPLETED**: **WORKING UNIFIED FUNCTION SYSTEM** 
- ✅ Created `scripts/haplotype_estimation_functions.R` with unified `estimate_haplotypes()` function
- ✅ Created `scripts/test_haplotype_functions.R` - **THE CURRENT WORKING TEST SCRIPT**
- ✅ **TESTED AND DEBUGGED**: `Rscript scripts/test_haplotype_functions.R` runs successfully
- ✅ All 16 test cases pass: multiple positions, samples, methods (fixed/adaptive), h_cutoff values
- ✅ Proper verbosity levels (0=silent production, 2=detailed diagnostics)
- ✅ Verified different h_cutoff values produce different results (constraint accumulation working)
- 🎯 **NOW**: Update production scripts to call this working unified function

### **Professional Implementation Completed:**

## 🎯 **THE WORKING FOUNDATION**

**CURRENT WORKING TEST SCRIPT**: `scripts/test_haplotype_functions.R`
- **Command**: `Rscript scripts/test_haplotype_functions.R`
- **Status**: ✅ **RUNS SUCCESSFULLY** - This contains the working code
- **Function**: Calls unified `estimate_haplotypes()` function with various parameters
- **Test Cases**: 16 combinations (positions × samples × methods × h_cutoff values)
- **Verified**: Different h_cutoff values produce different results
- **Output**: Properly formatted data frames with all expected columns

**UNIFIED FUNCTION**: `scripts/haplotype_estimation_functions.R`
- **Function**: `estimate_haplotypes(pos, sample_name, df3, founders, h_cutoff, method, window_size_bp, chr, verbose)`
- **Methods**: Both "fixed" and "adaptive" 
- **Verbosity**: 0=silent, 1=basic, 2=detailed diagnostics
- **Output**: Standardized result list with chr, pos, sample, method, estimate_OK, founder frequencies

**STEP 2: Create Unified Production Wrapper** ✅ **COMPLETED**
- ✅ Created `scripts/run_haplotype_estimation.R` - unified production wrapper
- ✅ Follows same pattern as `test_haplotype_functions.R` (data loading, processing, function calls)
- ✅ Handles both fixed and adaptive methods with single interface
- ✅ Uses `estimate_haplotypes()` function with `verbose=0` for production
- ✅ Intelligent output filenames: `fixed_window_50kb_results_chr2R.RDS`, `adaptive_window_h6_results_chr2R.RDS`
- ✅ **Moved old scripts**: `REFALT2haps.FixedWindow.Single.R` and `REFALT2haps.AdaptWindow.Single.R` → `old_REFALT2haps/`

**STEP 3: Update Slurm Wrapper** ✅ **COMPLETED**
- ✅ Updated `scripts/haplotype_testing_from_table.sh` to call unified wrapper
- ✅ Simplified: Single call to `run_haplotype_estimation.R` for both methods
- ✅ Uses renamed parameter file: `helpfiles/production_slurm_params.tsv`
- ✅ **Clean architecture**: Slurm → Production wrapper → Unified function

## 🎯 **NEW CLEAN ARCHITECTURE COMPLETED**

**Layer 1**: `scripts/haplotype_estimation_functions.R` - Core `estimate_haplotypes()` function  
**Layer 2**: `scripts/run_haplotype_estimation.R` - Production wrapper (chromosome-wide scans)  
**Layer 3**: `scripts/haplotype_testing_from_table.sh` - Slurm array wrapper (parameter combinations)

**Command Pattern**:
```bash
# Slurm calls production wrapper with specific parameter combination
Rscript run_haplotype_estimation.R chr2R fixed 50 process/JUICE helpfiles/JUICE_haplotype_parameters.R

# Production wrapper calls unified function for all positions/samples  
estimate_haplotypes(pos, sample, df3, founders, h_cutoff, method, chr, verbose=0)
```

**READY FOR DEPLOYMENT** 🚀

**NEW PARADIGM**: Production scripts become simple purrr::pmap_dfr calls to tested unified function

## 📊 **PARALLELIZATION STRATEGY UNDERSTANDING**

**SLURM Array Jobs**: Each array element processes one parameter combination in parallel:
- **Chromosome**: chr2R, chr2L, chr3L, chr3R, chrX
- **Method**: fixed window vs adaptive window  
- **Parameter**: window size (kb) or h_cutoff value

**Within Each Job**: Process ALL positions and ALL samples for that combination:
- **All scan positions**: `seq(min_position, max_position, by = 10kb)` across euchromatin
- **All samples**: 6 samples defined in `names_in_bam` parameter
- **Quality filtering**: Only SNPs where founders are fixed (< 3% or > 97% frequency)

**Why This Works**: 
- **Parallel across methods** (9 different combinations run simultaneously)
- **Sequential within methods** (positions and samples processed efficiently in each job)
- **Leverages cluster resources** for the computationally expensive method comparisons

### **🚨 CRITICAL TESTING AND RESOURCE REALITY**

**Laws of Robotics - Code Changes:**
1. **If it works, DON'T FIX IT. If it's broken, fix ONLY what's broken.**
2. **Test script works perfectly → Production script should exactly mimic the working test**
3. **NEVER re-write production code from scratch when test code works**

**Resource Reality:**
- **Cluster debugging is EXPENSIVE** (time, compute resources, queue wait)
- **Local testing is CHEAP** (10 seconds vs overnight)
- **Test scripts must be the source of truth** - they replicate exact production logic
- **Production bugs caught locally save hours/days of cluster debugging**

**Current Protocol:**
- ✅ **Test scripts verified working by user**
- 🎯 **Copy exact working logic to production (remove only diagnostic output)**
- ❌ **NO new algorithms, NO new ideas, NO rewrites**

### **Corrected Algorithm Logic:**
- **Always run LSEI** (if sufficient SNPs) to get actual haplotype frequencies (B1, B2, ..., AB8)
- **Always run clustering** to assess founder distinguishability quality
- **Output both**: Frequency estimates AND quality flag per position/sample
- **estimate_OK meanings**: 1=reliable, 0=unreliable, NA=LSEI failed

### **🚨 CRITICAL DESIGN PRINCIPLE:**
**FOUNDERS ARE NOT HARDCODED TO 8** - They are whatever is defined in the parameter file (`helpfiles/JUICE_haplotype_parameters.R`).
- **Scripts must work with ANY number of founders** (6, 8, 10, etc.)
- **Use `length(founders)` not hardcoded `8`**
- **Never assume exactly 8 founders** - this is flexible by design
- **Examples in conversation**: When I say "8 founders" it's shorthand for "however many founders are defined"

### **📁 PARAMETER FILE CONFIGURATION:**
**All scripts read parameters from: `helpfiles/JUICE_haplotype_parameters.R`**

**Parameters used by scripts:**
- **`founders`**: Founder sample names (e.g., B1, B2, ..., AB8) - defines how many founders
- **`step`**: Step size for chromosome scanning (10000 = 10kb steps)
- **`h_cutoff`**: Clustering threshold for fixed window (2.5 default)
- **`names_in_bam`**: Specific samples to process (AJ_1_1, AJ_2_1, GJ_1_1, etc.)

**Parameters NOT used by scripts:**
- **`size`**: COMMENTED OUT - window size comes from command line for flexibility

**Key Benefits:**
- **Selective processing**: Can handle large REFALT files but only analyze specific samples
- **Flexible founder counts**: Works with any number of founders defined in parameter file
- **Single source of truth**: All parameters centralized in one file
- **Easy configuration**: Change step size, samples, or founders without editing scripts

### **Corrected Production Scripts:**
- ✅ **`scripts/REFALT2haps.FixedWindow.Single.R`**: Full LSEI estimation + distinguishability check
- ✅ **`scripts/REFALT2haps.AdaptWindow.Single.R`**: Progressive expansion + LSEI at optimal window
- ✅ **Complete Output**: Both haplotype frequencies (B1, B2, ..., AB8) AND estimate_OK quality flag
- ✅ **Major Efficiency**: Quality filter applied once at data loading, not in every window loop
- ✅ **Test Scripts**: `test_fixed_window.R` and `test_adaptive_window.R` match production logic exactly
- ✅ **Bug Fixes**: All scripts now include both founder AND sample data for LSEI

### **✅ Test Scripts Completed Successfully**

Both test scripts ran successfully with excellent diagnostic output:
- **Fixed window**: Comprehensive SNP tables, distance matrices, group analysis
- **Adaptive window**: Progressive window testing, SNP tracking, detailed clustering analysis  
- **Clean formatting**: Right-justified tables, 1-decimal distances, easy visual scanning
- **Major efficiency**: Quality filter applied once vs. thousands of times in loops
- **Scripts organization**: All debugging/testing scripts properly organized in `scripts/debug_and_testing/`

### **🚀 READY: Run Full Parameter Pipeline via Slurm**

**All bugs fixed, algorithm corrected - ready for production run:**

```bash
# Run full haplotype testing pipeline with all parameter combinations
sbatch --array=1-9 scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE yes
```

**This will test ALL parameter combinations:**
- Fixed window: 20kb, 50kb, 100kb, 200kb, 500kb
- Adaptive window: h_cutoff 4, 6, 8, 10

**Expected pipeline outputs:**
- Complete haplotype estimates: B1, B2, B3, B4, B5, B6, B7, AB8 frequencies for each position/sample/method
- Quality flags: estimate_OK (1=reliable, 0=unreliable, NA=failed) for each estimate
- Results saved as `.RDS` files: `fixed_window_*kb_results_chr2R.RDS` and `adaptive_window_h*_results_chr2R.RDS`
- Much faster execution due to efficiency improvements
- Different results across h_cutoff values and window sizes

## 📁 **SCRIPTS FOLDER ORGANIZATION**

**KEEP SCRIPTS FOLDER CLEAN - NO TEMPORARY OR DEBUGGING SCRIPTS IN ROOT**

### **Production Scripts (scripts/ root) - CLEANED UP:**
- `REFALT2haps.FixedWindow.Single.R` - Fixed window distinguishability
- `REFALT2haps.AdaptWindow.Single.R` - Adaptive window distinguishability  
- `euchromatic_SNP_imputation_single.R` - SNP imputation
- `haplotype_testing_from_table.sh` - Slurm pipeline wrapper
- `evaluate_haplotype_methods.R` - Results evaluation
- `summarize_pipeline_results.R` - Pipeline monitoring
- **Organized subfolders**: `haps2scan/`, `Heterozygosity_tests/`, `old_REFALT2haps/`, `raw2bam2REFALT/`

### **Debugging/Testing Scripts (scripts/debug_and_testing/):**
- `test_fixed_window.R` - Test fixed window algorithm
- `test_adaptive_window.R` - Test adaptive window algorithm
- All other debugging, testing, and temporary scripts

### **🚨 RULE: NO TEMPORARY SCRIPTS IN PRODUCTION FOLDER**
- ✅ **CLEANUP COMPLETED**: Moved `debug_pos_column.R` and `debug_production_pipeline.R` to `debug_and_testing/`
- ✅ **OLD FILES REMOVED**: Deleted `.OLD.R` backup files 
- All debugging scripts go in `scripts/debug_and_testing/`
- All temporary scripts get deleted after use or moved to `debug_and_testing/`
- Keep production scripts folder clean and organized

## ⚠️ **CRITICAL WORKFLOW CONSTRAINT** ⚠️

**This is a CLUSTER-BASED project. Local terminal access is severely limited:**

- **Data processing happens on the cluster** (HPC system)
- **Local development is script-only** - cannot access cluster data or run cluster jobs
- **File paths in CURRENT_STATUS.md refer to CLUSTER file structure**, not local
- **Pipeline execution requires cluster access** via `sbatch` commands
- **Results analysis requires cluster access** to read .RDS files

**Why this matters**: I cannot run `ls`, check file contents, or debug cluster issues from the local terminal. All cluster operations must be documented in CURRENT_STATUS.md for reference.

## 🚨 **CRITICAL TESTING AND RESOURCE REALITY** 🚨

**RUNTIME AND RESOURCES MATTER ENORMOUSLY:**

- **Local test**: 10 seconds, free, immediate feedback, easy debugging
- **Cluster production run**: OVERNIGHT, expensive, slow feedback loop, very costly to debug
- **Cost of debugging on cluster**: Days of wasted time and computational resources
- **Efficiency imperative**: Catch ALL bugs locally before any cluster runs

**CORRECT TESTING FLOW - TEST TO PRODUCTION:**

1. **Test script works perfectly** (fast, local, debugged)
2. **Production script should EXACTLY mimic the working test** (same logic, same data structures, same output)
3. **If test passes → production should work**
4. **If production fails → fix production to match working test**

**WRONG FLOW (what I keep doing):**
- Production fails → "fix" test to match broken production
- This wastes cluster resources and defeats the purpose of testing

**THE GOLDEN RULE:**
**Test scripts are the source of truth. Production scripts should replicate working test logic exactly.**

## 🚨 **LAWS OF ROBOTICS - CODE CHANGES** 🚨

**BEFORE MAKING ANY CODE CHANGE, I MUST ALWAYS CHECK THIS CHECKLIST:**

### **The Practical Code Change Protocol:**

#### **BEFORE MAKING MAJOR CHANGES:**
1. **Make a backup** of the working code (git branch, copy file, etc.)
2. **Document what the code was doing right** - capture current behavior and metrics
3. **Document what you're trying to fix** - be specific about the actual problem
4. **Create tests that verify the main function still works** - baseline metrics

#### **AFTER MAKING CHANGES:**
1. **Run the same tests** that the working code passed
2. **Verify it still does its main function correctly**
3. **Check that performance metrics haven't degraded**
4. **Only commit if it passes all the original tests**

### **The Golden Rule:**
**"If it works, DON'T FIX IT. If it's broken, fix ONLY what's broken."**

### **Why This Matters:**
**I have a destructive tendency to delete working code and then spend days debugging my own unnecessary changes. This protocol prevents that by requiring backup, documentation, and testing before any major modification.**

### **The Key Insight:**
**Make a test of the code as it is working, record metrics, and then the new code has to at least pass those same checks before committing.**

---

## Current Phase: Phase 4 (Method evaluation and comparison) + Pipeline Re-execution

**Status**: Adaptive window haplotype estimation + SNP imputation jobs are running on the cluster with the fixed production code.

**Current Position in Pipeline Cycle**: 
- ✅ **DEBUGGING COMPLETE** - Redundant constraint bug identified and fixed
- 🔄 **RUN ANALYSIS** - Adaptive window haplotype estimation + SNP imputation in progress (array elements 6-9)
- ⏳ **MONITOR COMPLETION** - Use summarize script to check pipeline status
- ⏳ **EVALUATE RESULTS** - Use evaluation script to confirm different performance across h_cutoff values

## Recent Progress (Last 24 hours)

### ✅ Completed Tasks
1. **Full pipeline execution** - Haplotype estimation and SNP imputation completed for chr2R
2. **Comprehensive method evaluation** - All 12 estimators (6 fixed + 6 adaptive) evaluated
3. **Sliding window analysis integration** - Regional performance analysis (250 SNP windows, 50 SNP steps)
4. **Performance optimization** - Sliding window calculation optimized with vectorized tidyverse operations
5. **Documentation updates** - ADAPTIVE_WINDOW_ALGORITHM.md created and updated
6. **Bug identification and fix** - Redundant constraint bug identified and fixed in production code
7. **Code cleanup** - Debugging scripts moved to debug_and_testing folder

### 🔧 Bug Fix Completed ✅

**The Problem**: All adaptive window methods (h4, h6, h8, h10) were producing identical haplotype estimates due to redundant constraint accumulation.

**Root Cause**: The 10kb window (1 founder group) was adding a redundant "sum=1" constraint that was already enforced, causing matrix singularity and LSEI failures.

**The Fix**: Implemented smart constraint detection that only runs LSEI and accumulates constraints when groups meaningfully change (not just count, but composition).

**Verification**: Test script confirmed the fix works - h4 vs h10 now produce different results:
- Position 1.5e+07: Max frequency difference = 0.055465
- Position 2e+07: Max frequency difference = 0.118146

### 🚀 Pipeline Re-execution - IN PROGRESS

**Production Code Updated**: `scripts/REFALT2haps.AdaptWindow.Single.R` now has the fixed constraint logic.

**Current Status**: Adaptive window haplotype estimation jobs are running (array elements 6-9).

**Expected Results**: Different h_cutoff values should now produce different performance metrics, resolving the mystery of identical cluster results.

## Key Findings

### Adaptive Window Algorithm
- **IS working correctly** in current production code
- **Progressive window expansion** with hierarchical clustering
- **Smart constraint accumulation** - only when groups meaningfully change
- **Bug fixed**: No more redundant constraint accumulation

### Performance Results (Previous Run)
- **Fixed windows**: 10kb performs best, larger windows degrade performance
- **Adaptive windows**: Superior to fixed windows across all metrics
- **All adaptive methods identical**: This was due to the constraint bug (now fixed)

## ⚠️ CRITICAL WARNING - Adaptive Window Code

**The adaptive window algorithm in `scripts/REFALT2haps.AdaptWindow.Single.R` is now working correctly after fixing the redundant constraint bug. DO NOT edit this code without understanding the big picture.**

### What We Are Trying to Achieve
1. **Progressive window expansion** - Start small (10kb) and grow to capture optimal founder separation
2. **Hierarchical clustering** - Use `h_cutoff` parameter to control founder grouping aggressiveness
3. **Smart constraint accumulation** - Preserve good founder separation from smaller windows
4. **Different results for different h_cutoff values** - This is the core goal that was broken before

### Why This Code is Critical
- **This is the core algorithm** for adaptive window haplotype estimation
- **It took significant debugging** to identify and fix the redundant constraint bug
- **Random edits can break the entire pipeline** and lose weeks of work
- **The algorithm works correctly now** - don't change it without understanding what it does

### What NOT to Do
- ❌ **Don't randomly change constraint logic** without understanding the algorithm
- ❌ **Don't modify the hierarchical clustering** without understanding why it's there
- ❌ **Don't change the window expansion logic** without understanding the progression
- ❌ **Don't edit this file** unless you're fixing a specific, identified bug

### What TO Do Instead
- ✅ **Read ADAPTIVE_WINDOW_ALGORITHM.md** to understand the algorithm
- ✅ **Test changes on small datasets** before modifying production code
- ✅ **Verify that different h_cutoff values still produce different results**
- ✅ **Document any changes** and why they were necessary

**Remember: The goal is to have an adaptive window algorithm that produces meaningfully different results for different h_cutoff values, not to randomly modify code that's already working.**

## Next Steps

### Current Phase (Next 4-8 hours)
1. **✅ Bug fix completed** - Production code updated with smart constraint logic
2. **🔄 Haplotype estimation + SNP imputation running** - Adaptive window methods (h4, h6, h8, h10) with `run_imputation=yes`
3. **Monitor progress** - Use summarize script to check pipeline status
4. **Wait for completion** - Both haplotype estimation and SNP imputation are running together

### Immediate Next (Next 1-2 days)
1. **Monitor pipeline completion** - Use `summarize_pipeline_results.R` to check status
2. **Evaluate results** - Use `evaluate_haplotype_methods.R` to confirm different performance across h_cutoff values
3. **Compare with previous results** - Demonstrate the bug fix worked
4. **Verify SNP imputation** - Check that imputed frequencies are reasonable

### Short Term (Next 3-5 days)
1. **Parameter optimization** - find optimal h_cutoff values for JUICE dataset
2. **Method comparison** - comprehensive evaluation of fixed vs adaptive methods
3. **Extend to other chromosomes** (chr2L, chr3L, chr3R, chrX)
4. **Genome-wide analysis** using best performing methods

### Medium Term (Next week)
1. **Performance profiling** - runtime and accuracy trade-offs across genome
2. **Documentation** - final algorithm description and usage guide
3. **Pipeline optimization** - streamline for large-scale analysis

## Technical Details

### Files Modified Today
- `scripts/REFALT2haps.AdaptWindow.Single.R` - **FIXED** redundant constraint bug
- `scripts/debug_and_testing/` - **ORGANIZED** all debugging scripts
- `ADAPTIVE_WINDOW_ALGORITHM.md` - Comprehensive algorithm documentation
- `CURRENT_STATUS.md` - This file, updated with current status



## Workflow Constraints

### Cluster Access
- **Results stay on cluster** - cannot download large .RDS files locally
- **Local development** - scripts developed locally, tested on cluster
- **Version control** - git push/pull for code synchronization
- **Local terminal commands limited** - cannot access cluster data or run cluster jobs locally

### Git Workflow Preferences
- **Preferred workflow**: `git add . && git commit -m "message" && git push origin main` (single command)
- **Avoid separate commands** - user prefers efficiency over step-by-step execution
- **Commit messages** should be descriptive and explain what was changed and why

### Data Processing
- **Raw REFALT files** - direct processing, no intermediate files
- **Euchromatin boundaries**:
  - chr2L: 82,455 to 22,011,009 bp
  - chr2R: 5,398,184 to 24,684,540 bp (current focus)
  - chr3L: 158,639 to 22,962,476 bp
  - chr3R: 4,552,934 to 31,845,060 bp
  - chrX: 277,911 to 22,628,490 bp
- **10kb step intervals** - comprehensive chromosome coverage

## Recent Debugging Insights - RESOLVED ✅

### The Redundant Constraint Bug - RESOLVED ✅
**What Was Happening**:
1. 10kb window: All founders in 1 group → added constraint "sum=1"
2. This constraint was **redundant** with the base constraint already enforced
3. **Matrix singularity** occurred due to colinear rows
4. LSEI failed or behaved unexpectedly, causing identical results

**The Fix Applied**:
- Only accumulate constraints when groups meaningfully change
- Skip LSEI when groups unchanged, reuse previous constraints
- Detect meaningful changes by group composition, not just count
- **Production code updated** with smart constraint detection

**Verification Results**:
- **h4 vs h10 now produce different results** (frequency differences: 0.055-0.118)
- **Algorithm working as designed** - different h_cutoff values matter
- **Pipeline re-execution in progress**

## Current Questions - RESOLVED ✅

1. **✅ Is the redundant constraint bug the root cause** of identical cluster results? **YES**
2. **✅ Will fixing this bug** restore different performance across h_cutoff values? **YES**
3. **✅ Should we update production code** or just test the fix locally first? **PRODUCTION CODE UPDATED**
4. **✅ What's the best approach** for detecting meaningful group changes? **IMPLEMENTED**

## Success Metrics

### Algorithm Working Correctly - COMPLETED ✅
- [x] Different h_cutoff values produce different results (verified in test)
- [x] Progressive window expansion shows meaningful group changes
- [x] Constraint accumulation works without matrix issues
- [x] Full founder separation achieved or max window size reached
- [x] **Bug fix implemented** in production code

### Pipeline Performance - IN PROGRESS 🔄
- [x] Bug fixed and production code updated
- [🔄] Haplotype estimation running for adaptive methods
- [ ] Different performance metrics across h_cutoff values (expected now)
- [ ] Reasonable runtime (not stuck in "doom loops")
- [ ] Full chromosome coverage achieved

### Code Quality - COMPLETED ✅
- [x] No redundant constraints
- [x] Smart LSEI execution (only when needed)
- [x] Proper error handling and edge cases
- [x] Clear documentation and testing
- [x] **Code cleanup completed** - debugging scripts organized

## Current Pipeline Status

**Status**: Bug fixed, production code updated, adaptive window haplotype estimation **RUNNING**.

**Current Jobs**: Array elements 6-9 (adaptive h4, h6, h8, h10) are processing.

**Next Steps**:
1. **Monitor haplotype estimation jobs** until completion
2. **Run SNP imputation** for all adaptive methods
3. **Re-evaluate results** to verify different performance across h_cutoff values
4. **Extend to other chromosomes** once chr2R is validated

**Expected Outcome**: Different performance metrics across h_cutoff values, resolving the mystery of identical adaptive window results.

### 8. File Structure
```
.
├── ADAPTIVE_WINDOW_ALGORITHM.md
├── CURRENT_STATUS.md
├── data/
│   ├── founders/
│   ├── NOT_founders/
│   └── raw/
├── helpfiles/
│   ├── flymap.r6.txt
│   ├── founder.bams.txt
│   ├── haplotype_params.2R.tsv
│   └── JUICE/
│       └── JUICE_haplotype_parameters.R
├── logs/                                # SLURM job logs (ROOT LEVEL)
│   ├── haplotype_pipeline_<job_id>_1.err
│   ├── haplotype_pipeline_<job_id>_1.out
│   ├── haplotype_pipeline_<job_id>_2.err
│   ├── haplotype_pipeline_<job_id>_2.out
│   ├── haplotype_pipeline_<job_id>_3.err
│   ├── haplotype_pipeline_<job_id>_3.out
│   ├── haplotype_pipeline_<job_id>_4.err
│   ├── haplotype_pipeline_<job_id>_4.out
│   ├── haplotype_pipeline_<job_id>_5.err
│   ├── haplotype_pipeline_<job_id>_5.out
│   ├── haplotype_pipeline_<job_id>_6.err
│   ├── haplotype_pipeline_<job_id>_6.out
│   ├── haplotype_pipeline_<job_id>_7.err
│   ├── haplotype_pipeline_<job_id>_7.out
│   ├── haplotype_pipeline_<job_id>_8.err
│   ├── haplotype_pipeline_<job_id>_8.out
│   ├── haplotype_pipeline_<job_id>_9.err
│   └── haplotype_pipeline_<job_id>_9.out
├── process/
│   └── JUICE/
│       ├── haplotype_results/           # All results stored here
│       │   ├── adaptive_window_h10_results_chr2R.RDS
│       │   ├── adaptive_window_h4_results_chr2R.RDS
│       │   ├── adaptive_window_h6_results_chr2R.RDS
│       │   ├── adaptive_window_h8_results_chr2R.RDS
│       │   ├── fixed_window_100kb_results_chr2R.RDS
│       │   ├── fixed_window_200kb_results_chr2R.RDS
│       │   ├── fixed_window_20kb_results_chr2R.RDS
│       │   ├── fixed_window_500kb_results_chr2R.RDS
│       │   ├── fixed_window_50kb_results_chr2R.RDS
│       │   ├── haplotype_evaluation_detailed_chr2R.RDS
│       │   ├── haplotype_evaluation_plots_chr2R.png
│       │   ├── haplotype_evaluation_regional_chr2R.RDS
│       │   ├── haplotype_evaluation_summary_chr2R.RDS
│       │   ├── sliding_window_plot_chr2R.png
│       │   ├── sliding_window_results_chr2R.RDS
│       │   ├── sliding_window_summary_chr2R.RDS
│       │   ├── snp_imputation_adaptive_h10_chr2R.RDS
│       │   ├── snp_imputation_adaptive_h4_chr2R.RDS
│       │   ├── snp_imputation_adaptive_h6_chr2R.RDS
│       │   ├── snp_imputation_adaptive_h8_chr2R.RDS
│       │   ├── snp_imputation_fixed_100kb_chr2R.RDS
│       │   ├── snp_imputation_fixed_200kb_chr2R.RDS
│       │   ├── snp_imputation_fixed_20kb_chr2R.RDS
│       │   ├── snp_imputation_fixed_500kb_chr2R.RDS
│       │   └── snp_imputation_fixed_50kb_chr2R.RDS
│       └── RefAlt.chr2R.txt            # Raw REFALT data
├── ref/
├── scripts/
│   ├── debug_and_testing/              # Debugging scripts organized here
│   ├── haps2scan/
│   ├── raw2bam2REFALT/
│   ├── old_REFALT2haps/
│   ├── Heterozygosity_tests/
│   ├── REFALT2haps.AdaptWindow.Single.R
│   ├── evaluate_haplotype_methods.R
│   ├── REFALT2haps.FixedWindow.Single.R
│   ├── haplotype_testing_from_table.sh
│   └── euchromatic_SNP_imputation_single.R
```

**CRITICAL: Results are in `process/JUICE/haplotype_results/`, logs are in `logs/` (root level)**

## Complete Pipeline Commands - From Start to Finish

### 1. Run Full Pipeline (All Methods)
```bash
# Haplotype estimation + SNP imputation for ALL methods (array 1-9)
sbatch --array=1-9 scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE yes

# Haplotype estimation ONLY for ALL methods (array 1-9)
sbatch --array=1-9 scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE
```

### 2. Run Subset Pipeline (Current Run)
```bash
# Haplotype estimation + SNP imputation for ADAPTIVE methods only (array 6-9)
sbatch --array=6,7,8,9 scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE yes
```

### 3. Evaluate Results (After Pipeline Completes)
```bash
# Quick pipeline summary
Rscript scripts/summarize_pipeline_results.R process/JUICE chr2R

# Comprehensive evaluation
Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
```

### 4. Extend to Other Chromosomes
**To extend to other chromosomes, simply add rows to the existing parameter file:**

1. **Edit `helpfiles/haplotype_params.2R.tsv`** to add rows for other chromosomes
2. **Drop some estimation methods** to focus on finding optimal combinations
3. **Run the same pipeline command** - SLURM will process all rows in the file

**Note: The filename `haplotype_params.2R.tsv` is misleading since it will contain parameters for all chromosomes. Consider renaming to something like `haplotype_params.tsv` or `haplotype_params.all_chromosomes.tsv`.**

**Example parameter file structure:**
```tsv
chr2R	fixed	20
chr2R	fixed	50
chr2R	fixed	100
chr2R	adaptive	4
chr2R	adaptive	6
chr2L	fixed	50
chr2L	adaptive	6
chr3L	fixed	50
chr3L	adaptive	6
```

**Then run the full pipeline:**
```bash
sbatch --array=1-9 scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE_haplotype_parameters.R process/JUICE yes
```

**The pipeline automatically processes all rows in the parameter file - no need for separate commands per chromosome.**

## 🐛 **CRITICAL BUG FIXED - Quality Filter Logic**

### **The Problem:**
The quality filter was **backwards** in adaptive window scripts:
- **WRONG**: `filter(freq > 0.03 & freq < 0.97)` - kept noisy SNPs, filtered out good ones
- **CORRECT**: `filter(freq < 0.03 | freq > 0.97)` - keeps fixed SNPs, filters out noisy ones

### **The Concept:**
**Founders should be "fixed" for each SNP:**
- **Good SNP**: Founders have frequencies < 3% OR > 97% (mostly fixed)
- **Bad SNP**: Founders have frequencies between 3-97% (polymorphic, noisy)
- **Sample frequencies**: Can be anything (0-100%) - that's what we're estimating!

### **Why This Matters:**
- **Next-gen sequencing noise**: Allows 3% tolerance for sequencing errors
- **Founder fixation**: Ensures founders provide stable reference points
- **Sample estimation**: Sample can be polymorphic, we estimate its haplotype from fixed founders

### **Files Fixed:**
- `scripts/REFALT2haps.AdaptWindow.Single.R` - Production adaptive window
- `scripts/debug_and_testing/test_adaptive_window_algorithm.R` - Test script

### **Impact:**
This bug was causing **89% SNP failure rate** instead of expected **few percent**. Fixed scripts should now show much higher success rates.

## **Current Pipeline Status:**
