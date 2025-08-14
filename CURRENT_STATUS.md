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

## **CURRENT STATUS**
- ✅ Haplotype estimation scripts are now working correctly
- 🔄 Need to validate with small test case before running full pipeline
- 🎯 Big picture: Run haplotype testing for different parameter combinations on chr2R

## **NEXT STEPS**
1. **Quick validation test** (SMALL PICTURE):
   ```bash
   Rscript scripts/debug_haplotype_simple.R
   ```
   - Test 100 positions to verify scripts work correctly
   - Check output format and success rate
   - Should take ~1 minute, not hours

2. **If validation passes** (BIG PICTURE):
   ```bash
   sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE yes
   ```
   - Run full pipeline for all parameter combinations on chr2R
   - This is the actual goal

3. **Validate full results**:
   ```bash
   Rscript scripts/peek_haplotype_results.R process/JUICE/haplotype_results/fixed_window_20kb_results_chr2R.RDS
   ```

4. **Run evaluation script** to compare methods:
   ```bash
   Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
   ```

## **BIG PICTURE CONTEXT**
- **Project**: XQTL2 haplotype inference and SNP imputation
- **Dataset**: JUICE (Drosophila population)
- **Goal**: Compare fixed vs adaptive window methods for haplotype estimation
- **Ultimate objective**: Determine optimal parameters for genome-wide analysis

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
