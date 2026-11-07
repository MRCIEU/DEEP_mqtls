#!/usr/bin/env python3
"""DEEP global PCA v2 (01b_PCAplot): project a cohort onto the v2 reference PCA.

Called by 01b_PCAplot.sh. Projects a cleaned GRCh38 bfile with chr:pos_A1_A2 IDs onto
resources/global_pca/reference/deep_globalpca_v2_loadings.tsv.bgz and writes

  to --out-dir (uploaded, contains no sample IDs)
    <study>_globalPCA_v2_scores.tsv.gz    PC1-PC10 per sample, rows shuffled, no ID column
                                          (only when --share-level individual)
    <study>_globalPCA_v2_bins.tsv         sample counts in the fixed reference bin grid,
                                          counts of 1-4 reported as "<5"
    <study>_globalPCA_v2_coverage.tsv     reference SNPs found, allele swaps, mismatches, call rates
    <study>_globalPCA_v2_summary.tsv      sample size, coverage, PC means/SDs, versions
    <study>_globalPCA_v2.png              cohort over the reference (density only for aggregate)

  to --work-dir (kept by the cohort, never uploaded)
    <study>_globalPCA_v2_scores_with_ids.tsv.gz   FID, IID, PC1-PC20
    projected_ids.txt                              FID IID, used by the upload ID check

Alleles are matched to the reference by position; SNPs whose alleles are swapped are
recoded, SNPs whose alleles differ are dropped. Missing SNPs and missing calls contribute
zero (the reference mean), and genotypes are scaled by the full reference SNP count, as in
hail.experimental.pc_project, so scores are on the reference scale whatever the coverage.
"""

import argparse
import gzip
import os
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

PAIRS = [(1, 2), (1, 3), (2, 4)]
N_UPLOAD_PCS = 10
SUPPRESS_BELOW = 5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bfile", required=True)
    p.add_argument("--ref-dir", required=True)
    p.add_argument("--ref-prefix", default="deep_globalpca_v2")
    p.add_argument("--study", required=True)
    p.add_argument("--mode", required=True, choices=["post01b", "genetic_only"])
    p.add_argument("--share-level", required=True, choices=["individual", "aggregate"])
    p.add_argument("--build", default="38")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--min-coverage", type=float, default=0.95,
                   help="Stop if fewer than this fraction of reference SNPs are found (warn below 1.0)")
    p.add_argument("--pipeline-version", default="unknown")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mem-gb", type=int, default=16)
    return p.parse_args()


def project(args, ref):
    import hail as hl
    tmp = os.path.join(args.work_dir, "hail_tmp")
    os.makedirs(tmp, exist_ok=True)
    hl.init(quiet=True, local=f"local[{args.threads}]", log=os.path.join(args.work_dir, "hail_project.log"),
            tmp_dir=tmp, local_tmpdir=tmp, spark_conf={"spark.driver.memory": f"{args.mem_gb}g"})
    hl.default_reference("GRCh38")

    lo = hl.import_table(f"{ref}_loadings.tsv.bgz", delimiter="\t",
                         types={"locus": hl.tlocus("GRCh38"), "alleles": hl.tarray(hl.tstr),
                                "loadings": hl.tarray(hl.tfloat64), "pca_af": hl.tfloat64}).key_by("locus")
    lo = lo.checkpoint(os.path.join(tmp, "ref_loadings.ht"), overwrite=True)
    n_ref = lo.count()
    k = len(lo.aggregate(hl.agg.take(lo.loadings, 1))[0])

    mt = hl.import_plink(bed=f"{args.bfile}.bed", bim=f"{args.bfile}.bim", fam=f"{args.bfile}.fam",
                         reference_genome="GRCh38",
                         contig_recoding={str(i): f"chr{i}" for i in range(1, 23)}, skip_invalid_loci=True)
    n_samples = mt.count_cols()
    mt = mt.filter_rows(hl.is_defined(lo[mt.locus]))
    ra = lo[mt.locus].alleles
    mt = mt.annotate_rows(_direct=(mt.alleles[0] == ra[0]) & (mt.alleles[1] == ra[1]),
                          _swap=(mt.alleles[0] == ra[1]) & (mt.alleles[1] == ra[0]))
    rows = mt.rows()
    rstats = rows.aggregate(hl.struct(at_ref_loci=hl.agg.count(),
                                      direct=hl.agg.count_where(rows._direct),
                                      swapped=hl.agg.count_where(rows._swap)))
    mt = mt.filter_rows(mt._direct | mt._swap)
    # one variant per reference locus; drop loci matched more than once
    per_locus = mt.rows().group_by("locus").aggregate(n=hl.agg.count())
    mt = mt.filter_rows(per_locus[mt.locus].n == 1)
    mt = mt.annotate_entries(_GT=hl.if_else(mt._swap, hl.unphased_diploid_gt_index_call(2 - mt.GT.n_alt_alleles()),
                                            hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles())))
    mt = mt.annotate_rows(_cr=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt = mt.checkpoint(os.path.join(tmp, "cohort_matched.mt"), overwrite=True)
    cr = mt.aggregate_rows(hl.struct(n=hl.agg.count(), mean=hl.agg.mean(mt._cr), min=hl.agg.min(mt._cr),
                                     lt95=hl.agg.count_where(mt._cr < 0.95)))
    mt = mt.key_rows_by("locus")
    # loadings and af come from the reference table, so pc_project scales by all n_ref SNPs
    ht = hl.experimental.pc_project(mt._GT, lo.loadings, lo.pca_af)
    ht = ht.annotate(fid=mt.cols()[ht.s].fam_id)
    sc = ht.to_pandas()
    hl.stop()

    cov = {
        "reference_snps": n_ref,
        "cohort_variants_at_reference_loci": rstats.at_ref_loci,
        "matched_same_allele_order": rstats.direct,
        "matched_swapped_alleles": rstats.swapped,
        "allele_mismatch_dropped": rstats.at_ref_loci - rstats.direct - rstats.swapped,
        "matched_snps_used": cr.n,
        "reference_snps_absent": n_ref - cr.n,
        "coverage": cr.n / n_ref,
        "mean_call_rate_used_snps": cr.mean,
        "min_call_rate_used_snps": cr.min,
        "used_snps_call_rate_below_0.95": cr.lt95,
    }
    return sc, cov, n_samples, k


def bins_table(pcs, grid):
    out, clipped = [], {}
    for a, b in PAIRS:
        xa, xb = f"PC{a}", f"PC{b}"
        ga, gb = grid.loc[xa], grid.loc[xb]
        na = int(round((ga["max"] - ga["min"]) / ga.width))
        nb = int(round((gb["max"] - gb["min"]) / gb.width))
        ia = np.floor((pcs[xa] - ga["min"]) / ga.width).astype(int)
        ib = np.floor((pcs[xb] - gb["min"]) / gb.width).astype(int)
        clip = (ia < 0) | (ia >= na) | (ib < 0) | (ib >= nb)
        clipped[f"{xa}_{xb}"] = int(clip.sum())
        ia, ib = ia.clip(0, na - 1), ib.clip(0, nb - 1)
        for (i, j), n in pd.Series(1, index=pd.MultiIndex.from_arrays([ia, ib])).groupby(level=[0, 1]).sum().items():
            out.append({"x_pc": xa, "y_pc": xb,
                        "x_from": round(ga["min"] + i * ga.width, 4), "x_to": round(ga["min"] + (i + 1) * ga.width, 4),
                        "y_from": round(gb["min"] + j * gb.width, 4), "y_to": round(gb["min"] + (j + 1) * gb.width, 4),
                        "n": str(n) if n >= SUPPRESS_BELOW else f"<{SUPPRESS_BELOW}"})
    return pd.DataFrame(out), clipped


def plot(pcs, refd, study, share_level, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, INK2, MUTED, GRID, CTX, HI = "#1f1f1d", "#5f5e58", "#8a8983", "#e6e5e0", "#d6d5cf", "#2a78d6"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": GRID,
                         "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5})
    rows = ["points", "density"] if share_level == "individual" else ["density"]
    fig, axes = plt.subplots(len(rows), 3, figsize=(13, 4.3 * len(rows)), constrained_layout=True, squeeze=False)
    for c, (a, b) in enumerate(PAIRS):
        x, y = f"PC{a}", f"PC{b}"
        allx = pd.concat([refd[x], pcs[x]]); ally = pd.concat([refd[y], pcs[y]])
        rx, ry = allx.max() - allx.min(), ally.max() - ally.min()
        lim = ((allx.min() - rx * .05, allx.max() + rx * .05), (ally.min() - ry * .05, ally.max() + ry * .05))
        for r, kind in enumerate(rows):
            ax = axes[r, c]
            ax.scatter(refd[x], refd[y], s=3, c=CTX, linewidths=0, rasterized=True)
            for reg, s in refd.groupby("genetic_region"):
                ax.text(s[x].median(), s[y].median(), reg, fontsize=7, color=INK2, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.35))
            if kind == "points":
                ax.scatter(pcs[x], pcs[y], s=4, c=HI, linewidths=0, alpha=0.6, rasterized=True)
            else:
                cmap = matplotlib.colors.ListedColormap(plt.cm.Blues(np.linspace(0.4, 1.0, 256)))
                hb = ax.hexbin(pcs[x], pcs[y], gridsize=45, extent=(*lim[0], *lim[1]), mincnt=SUPPRESS_BELOW,
                               cmap=cmap, linewidths=0.2, edgecolors="white")
                if c == 2:
                    fig.colorbar(hb, ax=ax, shrink=0.8, label=f"samples per cell (cells < {SUPPRESS_BELOW} hidden)")
            ax.set_xlim(lim[0]); ax.set_ylim(lim[1])
            ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(length=2)
            ax.set_xlabel(x); ax.set_ylabel(y)
            if c == 0:
                ax.set_title({"points": "Samples", "density": "Density"}[kind], loc="left", fontsize=8.5, color=INK)
    fig.suptitle(f"{study}: projected onto the DEEP global PCA v2 reference (grey, HGDP + 1000 Genomes), "
                 f"n = {len(pcs):,}", x=0.01, ha="left", fontsize=10, color=INK)
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)
    ref = os.path.join(args.ref_dir, args.ref_prefix)
    stem = os.path.join(args.out_dir, f"{args.study}_globalPCA_v2")

    sc, cov, n_samples, k = project(args, ref)
    for key, val in cov.items():
        print(f"{key}: {val}", flush=True)
    if cov["coverage"] < args.min_coverage:
        sys.exit(f"Coverage {100 * cov['coverage']:.3f}% is below {100 * args.min_coverage:.0f}%: too few reference SNPs found. "
                 "Check the genome build and that SNP IDs are chr:pos_A1_A2. No outputs were written.")
    if cov["coverage"] < 1:
        print(f"WARNING: {cov['reference_snps_absent']} of {cov['reference_snps']} reference SNPs are absent "
              f"(coverage {100 * cov['coverage']:.3f}%); they contribute the reference mean", flush=True)

    S = np.vstack(sc.scores.values)
    pcs_all = pd.DataFrame(S, columns=[f"PC{i}" for i in range(1, k + 1)])
    local = pd.concat([sc[["fid", "s"]].rename(columns={"fid": "FID", "s": "IID"}).reset_index(drop=True),
                       pcs_all], axis=1)
    local.to_csv(os.path.join(args.work_dir, f"{args.study}_globalPCA_v2_scores_with_ids.tsv.gz"),
                 sep="\t", index=False, float_format="%.6f")
    local[["FID", "IID"]].to_csv(os.path.join(args.work_dir, "projected_ids.txt"), sep="\t", index=False, header=False)

    pcs = pcs_all[[f"PC{i}" for i in range(1, N_UPLOAD_PCS + 1)]]
    if args.share_level == "individual":
        shuffled = pcs.iloc[np.random.default_rng().permutation(len(pcs))].reset_index(drop=True)
        shuffled.to_csv(f"{stem}_scores.tsv.gz", sep="\t", index=False, float_format="%.6f")
    else:
        if os.path.exists(f"{stem}_scores.tsv.gz"):
            os.remove(f"{stem}_scores.tsv.gz")

    grid = pd.read_csv(f"{ref}_bins.tsv", sep="\t", comment="#").set_index("PC")
    bins, clipped = bins_table(pcs, grid)
    bins.to_csv(f"{stem}_bins.tsv", sep="\t", index=False)

    pd.DataFrame(list(cov.items()), columns=["metric", "value"]).to_csv(f"{stem}_coverage.tsv", sep="\t", index=False)

    summary = [("study", args.study), ("mode", args.mode), ("share_level", args.share_level),
               ("input_build", args.build), ("n_samples", n_samples), ("reference", args.ref_prefix),
               ("reference_snps", cov["reference_snps"]), ("coverage", round(cov["coverage"], 6)),
               ("pipeline_version", args.pipeline_version),
               ("completed_utc", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"))]
    summary += [(f"bins_clipped_{p}", n) for p, n in clipped.items()]
    for c in pcs.columns:
        summary += [(f"{c}_mean", round(pcs[c].mean(), 6)), (f"{c}_sd", round(pcs[c].std(), 6))]
    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(f"{stem}_summary.tsv", sep="\t", index=False)

    refd = pd.read_csv(f"{ref}_scores_flat.tsv.gz", sep="\t")
    plot(pcs, refd, args.study, args.share_level, f"{stem}.png")
    print(f"Wrote {stem}_*.tsv(.gz) and {stem}.png; IDs kept only in {args.work_dir}", flush=True)


if __name__ == "__main__":
    main()
