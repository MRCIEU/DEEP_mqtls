#!/usr/bin/env bash
# DEEP global PCA v2, step 2a: subset the phased HGDP+1kG reference to the
# reference samples and the step 1 common SNPs, and convert to one PLINK bfile.
#
# Usage: subset_reference.sh <ref_bcf_dir> <common.snplist.gz> <current_ref_scores.tsv> <out_dir> [threads]
#   ref samples = IDs in the current reference scores file that are present in the phased panel
set -euo pipefail

REF_BCF=$1; SNPLIST=$2; SCORES=$3; OUT=$4; THREADS=${5:-8}
BCFTOOLS=${BCFTOOLS:-bcftools}
PLINK2=${PLINK2:-resources/bin/plink2}
CHROMS=${CHROMS:-$(seq 1 22)}         # override for testing, e.g. CHROMS=22
bcf() { echo "$REF_BCF/hgdp1kgp_chr$1.filtered.SNV_INDEL.phased.shapeit5.bcf"; }

mkdir -p "$OUT/tmp"
tmp="$OUT/tmp"

# reference samples
first=$(echo $CHROMS | awk '{print $1}')
$BCFTOOLS query -l "$(bcf "$first")" | sort -u > "$tmp/panel_samples.txt"
tail -n +2 "$SCORES" | cut -f1 | sort -u > "$tmp/current_ref_samples.txt"
comm -12 "$tmp/current_ref_samples.txt" "$tmp/panel_samples.txt" > "$OUT/ref_samples.txt"
comm -23 "$tmp/current_ref_samples.txt" "$tmp/panel_samples.txt" > "$OUT/ref_samples_not_in_panel.txt"
echo "Reference samples: $(wc -l < "$OUT/ref_samples.txt") used, $(wc -l < "$OUT/ref_samples_not_in_panel.txt") not in panel"

# target sites (positions; alleles are matched later in make_reference_pca.py)
zcat "$SNPLIST" | awk -F'[:_]' -v OFS='\t' '{print "chr"$1, $2}' | sort -u -k1,1V -k2,2n > "$tmp/sites.tsv"
echo "Target sites: $(wc -l < "$tmp/sites.tsv")"

parts=()
for c in $CHROMS; do
    awk -v c="chr$c" '$1==c' "$tmp/sites.tsv" > "$tmp/sites_chr$c.tsv"
    [[ -s $tmp/sites_chr$c.tsv ]] || { echo "chr$c: no target sites, skipped"; continue; }
    $BCFTOOLS view --threads "$THREADS" -S "$OUT/ref_samples.txt" -T "$tmp/sites_chr$c.tsv" \
        -m2 -M2 -v snps -Ob -o "$tmp/chr$c.bcf" "$(bcf "$c")"
    $BCFTOOLS index -f "$tmp/chr$c.bcf"
    echo "chr$c: $($BCFTOOLS index -n "$tmp/chr$c.bcf") records"
    parts+=("$tmp/chr$c.bcf")
done

$BCFTOOLS concat --threads "$THREADS" -Ob -o "$tmp/ref_subset.bcf" "${parts[@]}"
$PLINK2 --bcf "$tmp/ref_subset.bcf" \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 \
    --output-chr 26 --make-bed --threads "$THREADS" \
    --out "$OUT/ref_subset"
echo "ref_subset: $(wc -l < "$OUT/ref_subset.bim") variants, $(wc -l < "$OUT/ref_subset.fam") samples"
rm -f "$tmp"/chr*.bcf "$tmp"/chr*.bcf.csi "$tmp"/sites_chr*.tsv
