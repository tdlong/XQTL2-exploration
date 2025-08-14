# CURRENT STATUS - XQTL2 Exploration Project

## **CURRENT GOAL**
Complete the haplotype estimation and SNP imputation pipeline for the JUICE dataset, with proper evaluation and validation of different methods.

## **CURRENT PHASE**
**Phase 3: Pipeline Validation and Optimization**
- ✅ Phase 1: Code organization and cleanup (COMPLETED)
- ✅ Phase 2: Bug fixes and debugging (COMPLETED)
- 🔄 Phase 3: Pipeline validation and optimization (IN PROGRESS)

## **WHAT WE JUST FIXED**
**Founder Column Count Bug** (Latest debugging session)
- **Problem**: All haplotype estimates were coming out as NA (0% success rate)
- **Root Cause**: Incorrect column count check in haplotype estimation scripts
  - Scripts were checking `ncol(founder_data) < length(founders) + 2`
  - Should be `ncol(founder_data) < length(founders) + 1` (for POS column)
- **Files Fixed**: 
  - `scripts/REFALT2haps.FixedWindow.Single.R`
  - `scripts/REFALT2haps.AdaptWindow.Single.R`
- **Validation**: Debug script confirmed 100% success rate after fix

## **CURRENT STATUS**
- ✅ Haplotype estimation scripts are now working correctly
- ✅ Pipeline is ready to run with proper parameters
- 🔄 Ready to test the full pipeline on chr2R

## **NEXT STEPS**
1. **Run the full pipeline** for chr2R:
   ```bash
   sbatch scripts/haplotype_testing_from_table.sh helpfiles/haplotype_params.2R.tsv helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE yes
   ```

2. **Validate results** using the peek script:
   ```bash
   Rscript scripts/peek_haplotype_results.R process/JUICE/haplotype_results/fixed_window_20kb_results_chr2R.RDS
   ```

3. **Run evaluation script** to compare methods:
   ```bash
   Rscript scripts/evaluate_haplotype_methods.R chr2R helpfiles/JUICE/JUICE_haplotype_parameters.R process/JUICE
   ```

4. **Extend to other chromosomes** once chr2R is working

## **BIG PICTURE CONTEXT**
- **Project**: XQTL2 haplotype inference and SNP imputation
- **Dataset**: JUICE (Drosophila population)
- **Goal**: Compare fixed vs adaptive window methods for haplotype estimation
- **Ultimate objective**: Determine optimal parameters for genome-wide analysis

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
