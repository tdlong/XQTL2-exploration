# CURRENT STATUS - XQTL2 Exploration Project

## **CURRENT GOAL**
Complete the haplotype estimation and SNP imputation pipeline for the JUICE dataset, with proper evaluation and validation of different methods.

## **CURRENT PHASE**
**Phase 3: Pipeline Validation and Optimization**
- ✅ Phase 1: Code organization and cleanup (COMPLETED)
- ✅ Phase 2: Bug fixes and debugging (COMPLETED)
- 🔄 Phase 3: Pipeline validation and optimization (IN PROGRESS)

## **WHAT WE JUST FIXED**
**Row Count Bug** (Latest debugging session)
- **Problem**: Success rate showing 800% due to incorrect row counting
- **Root Cause**: `founder_frequencies` list column causing `bind_rows()` to expand into multiple rows
- **Files Fixed**: 
  - `scripts/REFALT2haps.FixedWindow.Single.R`
  - `scripts/REFALT2haps.AdaptWindow.Single.R`
- **Fix**: Removed `founder_frequencies` list column, kept individual founder columns

**Founder Frequencies Column Artifact** (Latest fix)
- **Problem**: `founder_frequencies` column appearing as NA in output files
- **Root Cause**: Leftover references to `founder_frequencies` in result rows
- **Fix**: Removed all remaining `founder_frequencies` references from scripts

## **CURRENT STATUS**
- ✅ Haplotype estimation scripts are now working correctly
- ✅ Small test case validated (100% success rate, correct row counts)
- ✅ Full pipeline running successfully on chr2R
- ✅ **Fixed window haplotype estimation: COMPLETE** (all 5 files: 20kb, 50kb, 100kb, 200kb, 500kb)
- 🔄 **Adaptive window haplotype estimation: RUNNING** (no files yet - slower than fixed)
- 🔄 **SNP imputation: RUNNING** (no files yet - slower than haplotype estimation)
- 🎯 Pipeline progressing as expected (adaptive slower than fixed, SNP imputation slower than haplotype estimation)

## **NEXT STEPS**
1. **Monitor pipeline progress** (on cluster):
   ```bash
   # Check SLURM job status
   squeue
   
   # Check logs
   ls -la logs/
   
   # Check result files (from cluster tree output):
   # Fixed windows: ✓ All 5 complete (20kb, 50kb, 100kb, 200kb, 500kb)
   # Adaptive windows: ❌ Not yet complete
   # SNP imputation: ❌ Not yet complete
   ```
   
2. **Analysis preparation** (on cluster when ready):
   ```bash
   # Validate haplotype results
   Rscript scripts/peek_haplotype_results.R process/JUICE/haplotype_results/fixed_window_20kb_results_chr2R.RDS
   
   # Run method evaluation (when adaptive windows complete)
   Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
   ```
   
   - Haplotype estimation: Fixed windows complete, adaptive windows running
   - SNP imputation: Running for smaller fixed windows
   - Expected: Adaptive slower than fixed, SNP imputation slower than haplotype estimation

2. **Once pipeline completes**:
   ```bash
   Rscript scripts/peek_haplotype_results.R process/JUICE/haplotype_results/adaptive_window_h4_results_chr2R.RDS
   ```
   - Validate adaptive window results

3. **Run evaluation script** to compare methods:
   ```bash
   Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
   ```
   - Compare fixed vs adaptive window performance
   - Analyze SNP imputation accuracy

4. **Extend to other chromosomes** once chr2R analysis is complete

## **MONITORING & DEBUGGING STRATEGIES**
- **SLURM Job Monitoring**: Check job status with `squeue` and logs in `logs/` directory
- **Cluster File Monitoring**: Use `tree` or `ls` commands on cluster to check result files
- **Small-Scale Testing**: Use `test_haplotype_100.R` for quick validation on subset data
- **Result Validation**: Use `peek_haplotype_results.R` to inspect output files (on cluster)
- **Performance Profiling**: Use scripts in `scripts/debug_and_testing/` for detailed analysis
- **Incremental Development**: Test changes on small datasets before full pipeline runs
- **Analysis Commands**: All analysis happens on cluster - results never pulled locally

## **BIG PICTURE CONTEXT**
- **Project**: XQTL2 haplotype inference and SNP imputation
- **Dataset**: JUICE (Drosophila population)
- **Goal**: Compare fixed vs adaptive window methods for haplotype estimation
- **Ultimate objective**: Determine optimal parameters for genome-wide analysis

## **WORKFLOW CONSTRAINTS** ⚠️
- **SLURM Cluster Workflow**: Using high-performance cluster for computation, Cursor locally for development
- **Development Cycle**: git add + commit + push → pull on cluster → run SLURM jobs
- **Large Datasets**: Many steps take hours/days to complete on full datasets
- **Debugging Strategy**: Write small test scripts (e.g., `test_haplotype_100`) for subset testing
- **SLURM Scripts**: Main pipeline uses `scripts/haplotype_testing_from_table.sh` for job submission
- **Results Stay on Cluster**: Results are too large to pull back locally - all analysis happens on cluster
- **Local Access**: Only code, parameters, and small test files are available locally
- **Script Testing**: Any new R scripts require git add + commit + push → pull on cluster to test
- **Code Cleanup**: Deprecated scripts should be moved to `scripts/debug_and_testing/` or deleted

## **FUNDAMENTAL RULES** ⚠️
- **NEVER SKIP POSITIONS**: Every position/sample combination must have a result
- **Haplotype estimator returns**: Either founder frequency estimates OR NAs
- **No gaps in data**: If estimation fails, return NA values, don't skip
- **Complete coverage**: All positions × all samples must be represented in output

## **PIPELINE COMPONENTS**
1. **Haplotype Estimation**: Fixed window vs Adaptive window methods
2. **SNP Imputation**: Interpolate SNP frequencies from haplotype estimates
3. **Method Evaluation**: Compare accuracy and coverage of different approaches
4. **Validation**: Ensure results are biologically meaningful

## **KEY FILES**
- **Main Pipeline**: `scripts/haplotype_testing_from_table.sh`
- **Fixed Window**: `scripts/REFALT2haps.FixedWindow.Single.R`
- **Adaptive Window**: `scripts/REFALT2haps.AdaptWindow.Single.R`
- **SNP Imputation**: `scripts/euchromatic_SNP_imputation_single.R`
- **Evaluation**: `scripts/evaluate_haplotype_methods.R`
- **Debugging**: `scripts/peek_haplotype_results.R`

## **SUCCESS CRITERIA**
- [ ] Haplotype estimation produces >90% success rate
- [ ] SNP imputation produces reasonable frequency estimates
- [ ] Method evaluation shows meaningful differences between approaches
- [ ] Results are ready for genome-wide analysis

---
*Last Updated: After fixing founder column count bug*
*Next Update: After running full pipeline*
