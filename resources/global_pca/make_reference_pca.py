#!/usr/bin/env python
"""
DEEP global PCA v2, step 2b: build the reference PCA on SNPs shared by all cohorts.

Run once by the analysis team (not by cohorts). Takes HGDP+1kG reference
genotypes, restricts them to the reference samples and the step 1 common SNP
list, LD prunes, runs PCA, fixes the sign of each PC, and writes reference
loadings, allele frequencies and scores in the layout that
resources/datacheck/ancestry_infer.py already reads.

Defaults follow the current reference (hgdp_tgp_unrel_pass_filtGBMI_strictpruned):
MAF >= 0.001 in the reference samples, LD pruning r2 < 0.1 in 10 Mb windows,
long-range LD regions removed, 20 PCs from hl.hwe_normalized_pca.
Strand-ambiguous (A/T, C/G) SNPs are dropped unless --keep-ambiguous is set.

PC orientation: PC1 and PC2 are flipped so the EUR mean is negative (the rule
ancestry_infer.py applies at plot time); PC3 onwards are flipped so the region
with the largest absolute mean is positive. The same signs are applied to the
loadings and the scores, so projections need no further flipping.

Example:
  python resources/global_pca/make_reference_pca.py \
      --bfile step2/ref_subset \
      --snplist step1/deep_pca_v2.all_cohorts_common.snplist.gz \
      --samples step2/ref_samples.txt \
      --exclude-regions resources/genetics/exclude_b38.txt \
      --out-prefix step2/deep_globalpca_v2 \
      --threads 32 --mem-gb 200 --holdout-frac 0.1

Input genotypes: give exactly one of --mt, --vcf or --bfile (all GRCh38).
SNP list: one ID per line, chr:pos_A1_A2 with A1/A2 in ASCII order (pipeline format).
Samples: first column holds reference sample IDs (a header line is skipped).
"""

import argparse
import gzip
import json
import os
import random
import sys
from pathlib import Path

import hail as hl
import numpy as np
import pandas as pd

DEFAULT_META = (Path(__file__).resolve().parents[1] / "genetics" /
                "release_3.1.2_vcf_genomes_gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--mt", help="Hail MatrixTable of reference genotypes")
    src.add_argument("--vcf", help="Reference VCF/BCF-derived VCF(s), bgzipped; glob allowed")
    src.add_argument("--bfile", help="Reference PLINK bfile prefix (from subset_reference.sh)")
    p.add_argument("--snplist", required=True, help="Common SNP list (chr:pos_A1_A2)")
    p.add_argument("--samples", required=True, help="Reference sample IDs (first column)")
    p.add_argument("--exclude-regions", required=True, help="Long-range LD regions: chr start end [label]")
    p.add_argument("--sample-meta", default=str(DEFAULT_META), help="gnomAD HGDP+1kG sample metadata (for PC orientation)")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--n-pcs", type=int, default=20)
    p.add_argument("--maf", type=float, default=0.001, help="MAF filter in the reference samples")
    p.add_argument("--call-rate", type=float, default=0.99, help="Variant call rate filter in the reference samples")
    p.add_argument("--ld-r2", type=float, default=0.1)
    p.add_argument("--ld-window-bp", type=int, default=10_000_000)
    p.add_argument("--keep-ambiguous", action="store_true", help="Keep A/T and C/G SNPs (dropped by default)")
    p.add_argument("--holdout-frac", type=float, default=0.0,
                   help="If > 0, also rebuild the PCA without this fraction of samples, project them, and report shrinkage")
    p.add_argument("--seed", type=int, default=20261106)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--mem-gb", type=int, default=32)
    p.add_argument("--tmp-dir", default=None)
    return p.parse_args()


def read_ids(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        ids = [line.split()[0] for line in f if line.strip()]
    if ids and ids[0] in ("s", "IID", "ID", "sample", "#IID", "FID"):
        ids = ids[1:]
    return ids


def snp_key(contig, pos, alleles):
    """Pipeline-style ID chr:pos_A1_A2 with alleles in ASCII order, contig without 'chr'."""
    a = hl.sorted(alleles)
    return hl.str(contig).replace("^chr", "") + ":" + hl.str(pos) + "_" + a[0] + "_" + a[1]


def region_table(meta_path, sample_ids):
    with gzip.open(meta_path, "rt") as f:
        meta = pd.read_csv(f, sep="\t", usecols=["s", "hgdp_tgp_meta"])
    meta = meta[meta.s.isin(set(sample_ids))].copy()

    def get(m, key):
        try:
            return (json.loads(m) if isinstance(m, str) else m).get(key)
        except (TypeError, ValueError, AttributeError):
            return None
    meta["genetic_region"] = meta.hgdp_tgp_meta.apply(lambda m: get(m, "genetic_region"))
    meta["population"] = meta.hgdp_tgp_meta.apply(lambda m: get(m, "population"))
    return meta.drop(columns="hgdp_tgp_meta")


def orientation(scores_pd, regions):
    """Sign per PC: EUR mean negative for PC1-2; largest |region mean| positive for PC3+."""
    df = scores_pd.merge(regions, on="s", how="left")
    S = np.vstack(df.scores.values)
    signs, notes = [], []
    for j in range(S.shape[1]):
        means = pd.Series(S[:, j]).groupby(df.genetic_region.values).mean()
        if j < 2:
            flip = means.get("EUR", 0.0) > 0
            notes.append(f"PC{j + 1}: EUR mean {means.get('EUR', float('nan')):.4f} -> {'flipped' if flip else 'kept'}")
        else:
            top = means.abs().idxmax()
            flip = means[top] < 0
            notes.append(f"PC{j + 1}: largest |mean| {top} {means[top]:.4f} -> {'flipped' if flip else 'kept'}")
        signs.append(-1.0 if flip else 1.0)
    return signs, notes


def main():
    args = parse_args()
    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(out_dir, exist_ok=True)
    tmp_dir = args.tmp_dir or os.path.join(out_dir, "hail_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    hl.init(
        backend="spark",
        local=f"local[{args.threads}]",
        log=f"{args.out_prefix}.hail.log",
        tmp_dir=tmp_dir,
        local_tmpdir=tmp_dir,
        spark_conf={"spark.driver.memory": f"{args.mem_gb}g"},
    )
    hl.default_reference("GRCh38")

    summary = []

    def note(msg):
        print(msg, flush=True)
        summary.append(msg)

    # ---------------------------------------------------------------- inputs
    if args.mt:
        mt = hl.read_matrix_table(args.mt)
    elif args.vcf:
        mt = hl.import_vcf(args.vcf, reference_genome="GRCh38", force_bgz=True, array_elements_required=False)
    else:
        mt = hl.import_plink(bed=f"{args.bfile}.bed", bim=f"{args.bfile}.bim", fam=f"{args.bfile}.fam",
                             reference_genome="GRCh38",
                             contig_recoding={str(i): f"chr{i}" for i in range(1, 23)})
    mt = mt.select_rows().select_cols().select_entries("GT")

    sample_ids = read_ids(args.samples)
    note(f"Reference samples requested: {len(sample_ids)}")
    mt = mt.filter_cols(hl.literal(set(sample_ids)).contains(mt.s))

    snp_ids = read_ids(args.snplist)
    note(f"Common SNPs in input list: {len(snp_ids)}")
    snp_ht = hl.Table.parallelize([hl.struct(id=x) for x in snp_ids], hl.tstruct(id=hl.tstr),
                                  n_partitions=max(8, args.threads)).key_by("id")

    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
                        & mt.locus.in_autosome())
    mt = mt.annotate_rows(snp_id=snp_key(mt.locus.contig, mt.locus.position, mt.alleles))
    mt = mt.filter_rows(hl.is_defined(snp_ht[mt.snp_id]))
    mt = mt.checkpoint(os.path.join(tmp_dir, "ref_common.mt"), overwrite=True)
    n_var, n_samp = mt.count()
    note(f"Reference samples found: {n_samp}")
    note(f"Common SNPs found in reference: {n_var}")
    if n_samp < len(sample_ids):
        note(f"WARNING: {len(sample_ids) - n_samp} requested samples not found in genotype input")
    if n_var < len(snp_ids):
        note(f"WARNING: {len(snp_ids) - n_var} common SNPs not found in genotype input")

    if not args.keep_ambiguous:
        amb = hl.set([hl.set(["A", "T"]), hl.set(["C", "G"])])
        mt = mt.filter_rows(~amb.contains(hl.set(mt.alleles)))
        note(f"SNPs after dropping strand-ambiguous A/T and C/G: {mt.count_rows()}")

    intervals = []
    with open(args.exclude_regions) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3 or not parts[1].isdigit():
                continue
            c = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            intervals.append(hl.parse_locus_interval(f"{c}:{parts[1]}-{parts[2]}", reference_genome="GRCh38"))
    mt = hl.filter_intervals(mt, intervals, keep=False)
    note(f"SNPs after removing {len(intervals)} long-range LD regions: {mt.count_rows()}")

    # ------------------------------------------------------------- variant QC
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((mt.variant_qc.call_rate >= args.call_rate) & (hl.min(mt.variant_qc.AF) >= args.maf))
    mt = mt.drop("variant_qc")
    mt = mt.checkpoint(os.path.join(tmp_dir, "ref_qc.mt"), overwrite=True)
    note(f"SNPs after call rate >= {args.call_rate} and MAF >= {args.maf}: {mt.count_rows()}")

    # ------------------------------------------------------------ LD pruning
    pruned = hl.ld_prune(mt.GT, r2=args.ld_r2, bp_window_size=args.ld_window_bp)
    mt = mt.filter_rows(hl.is_defined(pruned[mt.row_key]))
    mt = mt.checkpoint(os.path.join(tmp_dir, "ref_pruned.mt"), overwrite=True)
    note(f"SNPs after LD pruning (r2 < {args.ld_r2}, window {args.ld_window_bp} bp): {mt.count_rows()}")

    # ------------------------------------------------------------------- PCA
    k = args.n_pcs
    mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=k, compute_loadings=True)
    loadings = loadings.annotate(pca_af=mt.rows()[loadings.key].pca_af)

    # fixed orientation, applied to loadings and scores together
    regions = region_table(args.sample_meta, sample_ids)
    scores_pd = scores.to_pandas()
    signs, sign_notes = orientation(scores_pd, regions)
    for n_ in sign_notes:
        note(f"Orientation {n_}")
    sl = hl.literal(signs)
    loadings = loadings.annotate(loadings=hl.range(k).map(lambda i: loadings.loadings[i] * sl[i]))
    scores = scores.annotate(scores=hl.range(k).map(lambda i: scores.scores[i] * sl[i]))

    # ---------------------------------------------------------------- export
    # same column layout as hgdp_tgp_unrel_pass_filtGBMI_strictpruned_*.tsv, block-gzipped
    loadings_out = f"{args.out_prefix}_loadings.tsv.bgz"
    scores_out = f"{args.out_prefix}_scores.tsv.bgz"
    loadings.select("loadings", "pca_af").export(loadings_out)
    scores.select("scores").export(scores_out)

    # flat table for plotting: s, region, population, PC1..PCk
    flat = scores.to_pandas()
    pcs = pd.DataFrame(np.vstack(flat.scores.values), columns=[f"PC{i}" for i in range(1, k + 1)])
    flat = pd.concat([flat[["s"]].reset_index(drop=True), pcs], axis=1).merge(regions, on="s", how="left")
    flat = flat[["s", "genetic_region", "population"] + [f"PC{i}" for i in range(1, k + 1)]]
    flat.to_csv(f"{args.out_prefix}_scores_flat.tsv.gz", sep="\t", index=False)

    snps = mt.rows().select("snp_id").to_pandas()
    snps["snp_id"].to_csv(f"{args.out_prefix}.snplist", index=False, header=False)

    with open(f"{args.out_prefix}_eigenvalues.txt", "w") as f:
        total = sum(eigenvalues)
        f.write("PC\teigenvalue\tprop_of_top_k\tsign\n")
        for i, ev in enumerate(eigenvalues, 1):
            f.write(f"PC{i}\t{ev}\t{ev / total}\t{int(signs[i - 1])}\n")

    # --------------------------------------------------- hold-out shrinkage
    if args.holdout_frac > 0:
        cols = sorted(mt.s.collect())
        held = set(random.Random(args.seed).sample(cols, max(1, int(len(cols) * args.holdout_frac))))
        hs = hl.literal(held)
        train = mt.filter_cols(~hs.contains(mt.s))
        test = mt.filter_cols(hs.contains(mt.s))
        train = train.annotate_rows(pca_af=hl.agg.mean(train.GT.n_alt_alleles()) / 2)
        _, tr_scores, tr_load = hl.hwe_normalized_pca(train.GT, k=k, compute_loadings=True)
        tr_load = tr_load.annotate(pca_af=train.rows()[tr_load.key].pca_af)
        proj = hl.experimental.pc_project(test.GT, tr_load.loadings, tr_load.pca_af)
        A = np.vstack(tr_scores.to_pandas().scores.values)
        B = np.vstack(proj.to_pandas().scores.values)
        ratio = B.std(axis=0) / A.std(axis=0)
        with open(f"{args.out_prefix}_holdout_shrinkage.tsv", "w") as f:
            f.write("PC\tsd_ratio_projected_over_in_sample\n")
            for i, r in enumerate(ratio, 1):
                f.write(f"PC{i}\t{r:.4f}\n")
        note(f"Hold-out ({len(held)} samples) SD ratio projected/in-sample, PC1-PC5: "
             + ", ".join(f"{r:.3f}" for r in ratio[:5]))

    note(f"Wrote {loadings_out}, {scores_out}, {args.out_prefix}_scores_flat.tsv.gz, {args.out_prefix}.snplist")
    with open(f"{args.out_prefix}_summary.txt", "w") as f:
        f.write("\n".join(summary) + "\n")
        f.write("Command: " + " ".join(sys.argv) + "\n")


if __name__ == "__main__":
    main()
