#!/usr/bin/env python3
"""DEEP global PCA v2, step 1d: summarise the common SNP list.

Reads the same cohort inputs as make_pca_common_snps.py (Haotian Tang) and
the reference SNP list it wrote, then reports for each cohort:
  - SNPs eligible and the fraction found in the reference (a low fraction
    points to a build or allele-coding problem)
  - SNPs that would be gained if that cohort were dropped (leave-one-out),
    which shows which cohorts limit the common list
It also checks that the common count matches Haotian's output and reports
counts per chromosome and the number of strand-ambiguous (A/T, C/G) SNPs.

Example
  python resources/global_pca/leave_one_out.py \
      --roots qc_links staging \
      --raw-ref step1/deep_pca_v2.raw.biallelic_snps.snplist.gz \
      --common step1/deep_pca_v2.all_cohorts_common.snplist.gz \
      --out-prefix step1/deep_pca_v2
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_pca_common_snps import cohort_inputs, id_key, read_cohort  # noqa: E402

AMBIGUOUS_CODES = {0 * 4 + 3, 1 * 4 + 2}  # A/T and C/G in Haotian's key encoding


def read_keys(path):
    with gzip.open(path, "rt") as f:
        return np.fromiter((id_key(line.strip()) for line in f), dtype=np.uint64)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--roots", nargs="+", required=True, help="Cohort roots, in the same order as given to make_pca_common_snps.py")
    p.add_argument("--raw-ref", required=True)
    p.add_argument("--common", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    ref = np.sort(read_keys(args.raw_ref))
    print(f"Reference SNPs: {len(ref)}", flush=True)
    cohorts, _ = cohort_inputs(args.roots)
    names, arrays, rows = [], [], []
    for name, path in sorted(cohorts.items()):
        keys, _ = read_cohort(path)
        arr = np.fromiter(keys, dtype=np.uint64, count=len(keys))
        del keys
        arr.sort()
        in_ref = arr[np.isin(arr, ref, assume_unique=True)]
        names.append(name)
        arrays.append(in_ref)
        rows.append({"cohort": name, "path": str(path), "eligible": len(arr), "in_reference": len(in_ref),
                     "frac_in_reference": len(in_ref) / len(arr) if len(arr) else 0.0})
        print(f"{name}: eligible {len(arr)}, in reference {len(in_ref)}", flush=True)

    n = len(arrays)
    u, cnt = np.unique(np.concatenate(arrays), return_counts=True)
    common = u[cnt == n]
    near = u[cnt == n - 1]
    for row, arr in zip(rows, arrays):
        row["gained_if_dropped"] = int(np.sum(~np.isin(near, arr, assume_unique=True)))

    haotian = read_keys(args.common)
    status = "match" if len(haotian) == len(common) and np.array_equal(np.sort(haotian), common) else "MISMATCH"
    print(f"Common SNPs: {len(common)} (make_pca_common_snps.py: {len(haotian)}, {status})", flush=True)

    with open(f"{args.out_prefix}.leave_one_out.tsv", "w") as out:
        cols = ["cohort", "eligible", "in_reference", "frac_in_reference", "gained_if_dropped", "path"]
        out.write("\t".join(cols) + "\n")
        for row in sorted(rows, key=lambda r: -r["gained_if_dropped"]):
            out.write("\t".join(f"{row[c]:.4f}" if c == "frac_in_reference" else str(row[c]) for c in cols) + "\n")

    chrom = (common >> np.uint64(4)) >> np.uint64(32)
    codes = common & np.uint64(15)
    ambiguous = int(np.isin(codes, list(AMBIGUOUS_CODES)).sum())
    with open(f"{args.out_prefix}.common_summary.tsv", "w") as out:
        out.write("item\tvalue\n")
        out.write(f"cohorts\t{n}\ncommon_snps\t{len(common)}\ncheck_vs_make_pca_common_snps\t{status}\n")
        out.write(f"strand_ambiguous_AT_CG\t{ambiguous}\ncommon_snps_non_ambiguous\t{len(common) - ambiguous}\n")
        for c in range(1, 23):
            out.write(f"chr{c}\t{int((chrom == c).sum())}\n")
    print(f"Strand-ambiguous: {ambiguous}; wrote {args.out_prefix}.leave_one_out.tsv and .common_summary.tsv", flush=True)
    if status != "match":
        sys.exit("Common list differs from make_pca_common_snps.py output; check the --roots order and inputs")


if __name__ == "__main__":
    main()
