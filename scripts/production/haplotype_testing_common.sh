#!/bin/bash

# Common functions for haplotype testing pipeline

function run_pipeline() {
    mkdir -p logs

    # Input Validation
    if [ ! -f "$PARAMS_TSV" ]; then
        echo "❌ Parameter TSV not found: $PARAMS_TSV"
        exit 1
    fi
    if [ ! -f "$PARFILE" ]; then
        echo "❌ Parameter file not found: $PARFILE"
        exit 1
    fi

    # Euchromatin Boundaries (Release 6 coordinates)
    declare -A euchromatin_start
    declare -A euchromatin_end

    euchromatin_start["chr2L"]=82455
    euchromatin_end["chr2L"]=22011009
    euchromatin_start["chr2R"]=5398184
    euchromatin_end["chr2R"]=24684540
    euchromatin_start["chr3L"]=158639
    euchromatin_end["chr3L"]=22962476
    euchromatin_start["chr3R"]=4552934
    euchromatin_end["chr3R"]=31845060
    euchromatin_start["chrX"]=277911
    euchromatin_end["chrX"]=22628490

    # Parse Parameters from TSV
    LINE_NUM=$SLURM_ARRAY_TASK_ID
    ROW=$(sed -n "${LINE_NUM}p" "$PARAMS_TSV")
    if [ -z "$ROW" ]; then
        echo "❌ No row found for array index $SLURM_ARRAY_TASK_ID"
        exit 1
    fi

    CHR=$(echo "$ROW" | awk -F"\t" '{print $1}')
    METHOD=$(echo "$ROW" | awk -F"\t" '{print $2}')
    PARAM=$(echo "$ROW" | awk -F"\t" '{print $3}')

    echo "=== Haplotype Pipeline Job ==="
    echo "Array ID: $SLURM_ARRAY_TASK_ID"
    echo "Chromosome: $CHR"
    echo "Method: $METHOD"
    echo "Parameter: $PARAM"
    echo "Parameter file: $PARFILE"
    echo "Output directory: $OUTDIR"
    echo "Run SNP imputation: $RUN_IMPUTATION"
    echo "Euchromatin region: ${euchromatin_start[$CHR]} - ${euchromatin_end[$CHR]} bp"
    echo "Time: $(date)"
    echo ""

    # Run haplotype estimation
    echo "=== Step 1: Haplotype Estimation ==="
    echo "Method: $METHOD"
    echo "Parameter: $PARAM"
    echo ""

    # Create results subdirectory
    RESULTS_DIR="$OUTDIR/haplotype_results"
    mkdir -p "$RESULTS_DIR"

    # Run haplotype estimation using unified production wrapper
    echo "Running haplotype estimation: $METHOD method with parameter $PARAM..."
    Rscript scripts/run_haplotype_estimation.R "$CHR" "$METHOD" "$PARAM" "$OUTDIR" "$PARFILE"
    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "✓ Haplotype estimation completed successfully"
    else
        echo "✗ Haplotype estimation failed"
        exit 1
    fi
}

