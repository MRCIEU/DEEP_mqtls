#!/usr/bin/env bash
# DEEP global PCA v2, step 2: reference PCA on the step 1 common SNPs.
#
# Run on ing-p6 from the repo root with deep_env active, after run_step1.sh:
#     PREFIX=deep_pca_v2_20261106 bash resources/global_pca/run_step2.sh
#
# 2a  subset the phased reference to the reference samples and common SNPs (subset_reference.sh)
# 2b  drop A/T and C/G SNPs, LD prune, PCA, orientation, hold-out test (make_reference_pca.py)
set -euo pipefail

ROOT=${ROOT:-/mnt/sod/DEEP/global_pca_v2}
REF_BCF=${REF_BCF:-$ROOT/ref_bcf}
PREFIX=${PREFIX:?set PREFIX to the step 1 prefix, e.g. deep_pca_v2_20261106}
COMMON=${COMMON:-$ROOT/step1/$PREFIX.all_cohorts_common.snplist.gz}
OUT=${OUT:-$ROOT/step2}
HAIL_TMP=${HAIL_TMP:-$OUT/hail_tmp}    # put on local disk if /mnt/sod is slow
PYTHON=${PYTHON:-python}
THREADS=${THREADS:-32}
MEM_GB=${MEM_GB:-200}
HOLDOUT=${HOLDOUT:-0.1}
EXTRA=${EXTRA:-}                       # extra make_reference_pca.py flags, e.g. "--keep-ambiguous"

here=resources/global_pca
[[ -f $here/make_reference_pca.py ]] || { echo "Run from the DEEP_mqtls repo root"; exit 1; }
[[ -s $COMMON ]] || { echo "Missing $COMMON (run step 1 first)"; exit 1; }
mkdir -p "$OUT"
exec > >(tee -a "$OUT/$PREFIX.run_step2.log") 2>&1
echo "== step 2 started $(date)  COMMON=$COMMON OUT=$OUT"

if [[ ! -s $OUT/ref_subset.bed ]]; then
    bash $here/subset_reference.sh "$REF_BCF" "$COMMON" \
        resources/genetics/hgdp_tgp_unrel_pass_filtGBMI_strictpruned_scores.tsv "$OUT" "$THREADS"
else
    echo "2a: $OUT/ref_subset.bed exists, skipping (delete it to redo)"
fi

$PYTHON $here/make_reference_pca.py \
    --bfile "$OUT/ref_subset" \
    --snplist "$COMMON" \
    --samples "$OUT/ref_samples.txt" \
    --exclude-regions resources/genetics/exclude_b38.txt \
    --out-prefix "$OUT/deep_globalpca_v2" \
    --threads "$THREADS" --mem-gb "$MEM_GB" --tmp-dir "$HAIL_TMP" \
    --holdout-frac "$HOLDOUT" $EXTRA

echo "== step 2 finished $(date)"
echo "Review: $OUT/deep_globalpca_v2_summary.txt, _eigenvalues.txt, _holdout_shrinkage.tsv"
