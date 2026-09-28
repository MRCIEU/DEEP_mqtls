#!/usr/bin/env bash
# Called by 01d/01g with parameters from config; GCTA validates its inputs.
# Records successful results in <out>.gwas_result.tsv. No plotting or GRM conversion.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash resources/genetics/run_gwas_with_fallback.sh \
  --gcta PATH --bfile PREFIX --pheno FILE --qcovar FILE --out PREFIX \
  --method auto|mlma|lr [--grm-sparse PREFIX] [--grm-dense PREFIX] \
  [--covar FILE] [--covar-maxlevel 500] [--thread-num 1]

auto: fastGWA-MLM first; MLMA only for the exact REML convergence error.
mlma: MLMA directly (for rerunning the other 01g groups after fallback).
lr:   fastGWA linear regression; never fall back to MLMA.

Use the matching thresholdeddense GRM for --grm-dense.
fastGWA uses --h2-limit 100; MLMA keeps the default 100 iterations.
Logs/results: <out>.run.*/. Manifest columns: method, result, run_dir.
EOF
}

die() { echo "Error: $*" >&2; exit 1; }

gcta=""
bfile=""
pheno=""
qcovar=""
out=""
method="auto"
sparse_grm=""
dense_grm=""
covar=""
maxlevel=500
threads=1

while (( $# )); do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --gcta) gcta=$2 ;;
        --bfile) bfile=$2 ;;
        --pheno) pheno=$2 ;;
        --qcovar) qcovar=$2 ;;
        --out) out=$2 ;;
        --method) method=$2 ;;
        --grm-sparse) sparse_grm=$2 ;;
        --grm-dense) dense_grm=$2 ;;
        --covar) covar=$2 ;;
        --covar-maxlevel) maxlevel=$2 ;;
        --thread-num) threads=$2 ;;
        *) die "Unknown option: $1" ;;
    esac
    shift 2
done

[[ -n "$out" ]] || die "--out is required"
[[ "$method" == auto || "$method" == mlma || "$method" == lr ]] || die "Invalid --method"

run_dir=$(mktemp -d "${out}.run.XXXXXX")
manifest="${out}.gwas_result.tsv"
if [[ -e "$manifest" ]]; then
    mv -- "$manifest" "$run_dir/previous.gwas_result.tsv"
fi
echo "GWAS run directory: $run_dir"

# Archive previously published results/plots before a new attempt starts.
# A failed rerun must not leave old outputs looking like a successful run.
for suffix in .fastGWA .fastGWA.gz .mlma .mlma.gz \
    _manhattan.pdf _nocisChr_manhattan.pdf _qqplot.jpeg _nocisChr_qqplot.jpeg _manhattan_beta.pdf; do
    if [[ -e "${out}${suffix}" ]]; then
        mv -- "${out}${suffix}" "$run_dir/previous${suffix}"
    fi
done

common=(--bfile "$bfile" --pheno "$pheno" --qcovar "$qcovar" --thread-num "$threads")
if [[ -n "$covar" ]]; then common+=(--covar "$covar"); fi

run_gcta() {
    local prefix=$1
    shift
    local -a statuses
    if "$gcta" "${common[@]}" "$@" --out "$prefix" 2>&1 | tee "${prefix}.console.log"; then
        return 0
    else
        statuses=("${PIPESTATUS[@]}")
        # A log-write error is not a REML failure and must not trigger fallback.
        (( statuses[1] == 0 )) || die "Could not write GCTA console log"
        return "${statuses[0]}"
    fi
}

validate_result() {
    local file=$1 expected=$2
    [[ -s "$file" ]] || die "Missing or empty result: $file"
    awk -v expected="$expected" '
        NR == 1 {
            count = split(expected, names, " ")
            if (NF != count) exit 1
            for (i = 1; i <= NF; i++) if ($i != names[i]) exit 1
            next
        }
        NF { if (NF != count) exit 1; rows++ }
        END { if (!rows) exit 1 }
    ' "$file" || die "Unexpected result header or missing/malformed rows: $file"
}

if [[ "$method" == auto || "$method" == lr ]]; then
    fast_args=(--fastGWA-lr)
    if [[ "$method" == auto ]]; then
        fast_args=(--fastGWA-mlm --grm-sparse "$sparse_grm")
    fi
    fast_prefix="$run_dir/fastgwa"
    if run_gcta "$fast_prefix" "${fast_args[@]}" --h2-limit 100 --covar-maxlevel "$maxlevel"; then
        final_method=fastGWA-lr
        if [[ "$method" == auto ]]; then final_method=fastGWA-mlm; fi
        result="${fast_prefix}.fastGWA"
        validate_result "$result" "CHR SNP POS A1 A2 N AF1 BETA SE P"
    else
        status=$?
        if [[ "$method" == auto ]] && grep -Fxq "Error: fastGWA-REML can't converge." "${fast_prefix}.console.log"; then
            echo "Detected exact fastGWA-REML convergence error; switching to MLMA."
            method=mlma
        else
            echo "GCTA failed (exit $status); no MLMA fallback. Logs: $run_dir" >&2
            exit "$status"
        fi
    fi
fi

if [[ "$method" == mlma ]]; then
    # Do not let MLMA compute a different GRM when the caller omitted its prefix.
    [[ -n "$dense_grm" ]] || die "MLMA requires --grm-dense"
    mlma_prefix="$run_dir/mlma"
    if run_gcta "$mlma_prefix" --mlma --grm "$dense_grm"; then
        final_method=mlma
        result="${mlma_prefix}.mlma"
        validate_result "$result" "Chr SNP bp A1 A2 Freq b se p"
    else
        status=$?
        echo "MLMA failed (exit $status); logs: $run_dir" >&2
        exit "$status"
    fi
fi

printf 'method\tresult\trun_dir\n%s\t%s\t%s\n' "$final_method" "$result" "$run_dir" > "$run_dir/result.tsv"
cp -- "$run_dir/result.tsv" "${manifest}"
echo "GWAS completed: method=$final_method result=$result"
