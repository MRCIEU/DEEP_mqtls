#!/usr/bin/env bash
# DEEP global PCA v2, step 1: SNP list common to all cohorts and the reference.
#
# Run on ing-p6 from the repo root with deep_env active:
#     bash resources/global_pca/run_step1.sh                 # all parts
#     STEPS="1b" bash resources/global_pca/run_step1.sh      # only convert the shared SNP lists
#     STEPS="1a 1c 1d" PREFIX=deep_pca_v2_run2 bash resources/global_pca/run_step1.sh
#
# 1a  link each QC-submitted cohort's results/01/data.afreq.gz into QC_LINKS
# 1b  convert the cohorts that shared SNP lists only (prepare_cohort_snplist.py)
# 1c  intersect cohorts and reference (make_pca_common_snps.py, Haotian Tang)
# 1d  leave-one-out and summary (leave_one_out.py)
#
# Every path can be overridden with an environment variable (see below).
set -euo pipefail

ROOT=${ROOT:-/mnt/sod/DEEP/global_pca_v2}
EXTRACTED=${EXTRACTED:-/mnt/sod/DEEP/sftp_downloads/extracted}
REF_BCF=${REF_BCF:-$ROOT/ref_bcf}
RAW=${RAW:-$ROOT/cohort_snplists}
STAGING=${STAGING:-$ROOT/staging}
QC_LINKS=${QC_LINKS:-$ROOT/qc_links}
OUT=${OUT:-$ROOT/step1}
PREFIX=${PREFIX:-deep_pca_v2_$(date +%Y%m%d_%H%M)}   # 1c never overwrites, so each run needs a new prefix
RAW_SNPLIST=${RAW_SNPLIST:-}          # reuse a previous <prefix>.raw.biallelic_snps.snplist.gz (skips the slow BCF scan)
STEPS=${STEPS:-1a 1b 1c 1d}
CONVERT=${CONVERT:-UGR QBB BMC_ACTREC CRELES Pelotas}
PYTHON=${PYTHON:-python}
BCFTOOLS=${BCFTOOLS:-bcftools}
THREADS=${THREADS:-8}
OVERWRITE=${OVERWRITE:-0}             # 1 = redo the 1b conversions

here=resources/global_pca
[[ -f $here/make_pca_common_snps.py ]] || { echo "Run from the DEEP_mqtls repo root"; exit 1; }
has() { [[ " $STEPS " == *" $1 "* ]]; }
mkdir -p "$OUT" "$STAGING"
log="$OUT/${PREFIX}.run_step1.log"
exec > >(tee -a "$log") 2>&1
echo "== step 1 ($STEPS) started $(date)"
echo "ROOT=$ROOT EXTRACTED=$EXTRACTED REF_BCF=$REF_BCF RAW=$RAW STAGING=$STAGING OUT=$OUT PREFIX=$PREFIX CONVERT=$CONVERT"

# ---- 1a. QC-submitted cohorts: newest results/01/data.afreq.gz per cohort folder
if has 1a; then
    echo "== 1a: linking QC-submitted cohorts from $EXTRACTED"
    rm -rf "$QC_LINKS"; mkdir -p "$QC_LINKS"
    links="$OUT/${PREFIX}.qc_links.tsv"
    printf "cohort\tstatus\tpath\n" > "$links"
    find -L "$EXTRACTED" -path '*/results/01/data.afreq.gz' -printf '%T@\t%p\n' | sort -rn |
    while IFS=$'\t' read -r _ path; do
        rel=${path#"$EXTRACTED"/}
        cohort=${rel%%/*}
        if [[ -e "$QC_LINKS/$cohort/results/01/data.afreq.gz" ]]; then
            printf "%s\tskipped_older_copy\t%s\n" "$cohort" "$path" >> "$links"
            continue
        fi
        mkdir -p "$QC_LINKS/$cohort/results/01"
        ln -s "$path" "$QC_LINKS/$cohort/results/01/data.afreq.gz"
        printf "%s\tused\t%s\n" "$cohort" "$path" >> "$links"
    done
    for d in "$EXTRACTED"/*/; do
        c=$(basename "$d")
        [[ -e "$QC_LINKS/$c" ]] || printf "%s\tno_data.afreq.gz\t%s\n" "$c" "$d" >> "$links"
    done
    echo "QC cohorts linked: $(grep -c $'\tused\t' "$links" || true); without data.afreq.gz: $(grep -c $'\tno_data' "$links" || true) (see $links)"
fi

# ---- 1b. cohorts that shared SNP lists only
if has 1b; then
    echo "== 1b: converting shared SNP lists from $RAW: $CONVERT"
    prep() {   # prep <cohort> <args...>
        local cohort=$1; shift
        if [[ -s "$STAGING/$cohort/results/01/data.afreq.gz" && $OVERWRITE != 1 ]]; then
            echo "$cohort: already converted, skipping (OVERWRITE=1 to redo)"; return
        fi
        local ow=()
        [[ $OVERWRITE == 1 ]] && ow=(--overwrite)
        $PYTHON $here/prepare_cohort_snplist.py --cohort "$cohort" --out-root "$STAGING" \
            --threads "$THREADS" "${ow[@]}" "$@"
    }
    for c in $CONVERT; do
        case $c in
            UGR)        prep UGR        --format ugr     --input "$RAW/UGR_DEEP_PCA.zip" --n-obs-max 1700 ;;
            QBB)        prep QBB        --format qbb     --input "$RAW/QBB_extracted_snp_freq.txt.zip" ;;
            BMC_ACTREC) prep BMC_ACTREC --format actrec  --input "$RAW/BMC_ACTREC.zip" --build 37 --n-obs-max 8682 \
                            --files "freq=FREQ.frq_s.csv,info=InfoScores_NoHLA_WithBasePairs-n-Chrs_s.csv" ;;
            CRELES)     prep CRELES     --format afreq   --input "$RAW/CRELES_data.afreq.gz" --n-obs-max 2130 ;;
            # Pelotas: corrected hct05 files (6 Nov), extracted from the encrypted zip into $RAW/Pelotas_DEEP512
            Pelotas)    prep Pelotas    --format pelotas --input "$RAW/Pelotas_DEEP512" --n-obs-max 1024 \
                            --files "freq=_hct05_freq.afreq,bim=_hct05.bim" ;;
            *) echo "No conversion rule for $c; add one here"; exit 1 ;;
        esac
    done
    # cohorts in STAGING but not in CONVERT would still be read by 1c: flag them
    for d in "$STAGING"/*/; do
        [[ -d $d ]] || continue
        c=$(basename "$d")
        [[ " $CONVERT " == *" $c "* ]] || echo "WARNING: $STAGING/$c is not in CONVERT but will be included in 1c; move it out to exclude it"
    done
fi

# ---- 1c. intersection (Haotian's script; changes are marked 'DEEP global PCA v2 change')
if has 1c; then
    [[ -d $QC_LINKS ]] || { echo "Missing $QC_LINKS: run 1a first"; exit 1; }
    for c in $(seq 1 22); do
        f="$REF_BCF/hgdp1kgp_chr${c}.filtered.SNV_INDEL.phased.shapeit5.bcf"
        [[ -s $f && -s $f.csi ]] || { echo "Missing $f or its .csi. Download into $REF_BCF with:"
            echo '  base=https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2'
            echo '  for c in $(seq 1 22); do f=hgdp1kgp_chr${c}.filtered.SNV_INDEL.phased.shapeit5.bcf; wget -c $base/$f $base/$f.csi; done'
            exit 1; }
    done
    echo "== 1c: intersecting cohorts with the reference"
    $PYTHON -B -u $here/make_pca_common_snps.py \
        --bcf-dir "$REF_BCF" \
        --pipeline-qc-dir "$QC_LINKS" \
        --pipeline-debug-dir "$STAGING" \
        --bcftools "$BCFTOOLS" \
        --output-dir "$OUT" \
        --prefix "$PREFIX" \
        ${RAW_SNPLIST:+--raw-snplist "$RAW_SNPLIST"}
fi

# ---- 1d. leave-one-out and summary
if has 1d; then
    echo "== 1d: leave-one-out"
    $PYTHON $here/leave_one_out.py \
        --roots "$QC_LINKS" "$STAGING" \
        --raw-ref "${RAW_SNPLIST:-$OUT/$PREFIX.raw.biallelic_snps.snplist.gz}" \
        --common "$OUT/$PREFIX.all_cohorts_common.snplist.gz" \
        --out-prefix "$OUT/$PREFIX"
    echo "Common SNP list: $OUT/$PREFIX.all_cohorts_common.snplist.gz"
    echo "Review: $OUT/$PREFIX.cohort_inputs.tsv, .leave_one_out.tsv, .common_summary.tsv"
fi
echo "== step 1 finished $(date)  PREFIX=$PREFIX"
