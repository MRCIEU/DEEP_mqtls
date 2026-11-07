#!/usr/bin/env python3
"""DEEP global PCA v2: check results/01b_pca before upload. Exits 1 on any problem.

1. Only expected files are present (anything else must be moved out before upload).
2. The scores file holds exactly PC1-PC10 as numbers: no ID or other column.
3. The bins file holds only the bin columns, with counts >= 5 or "<5".
4. No sample ID (FID or IID) appears as a whole token in any text file, including logs.
   IDs that are only digits and shorter than 6 characters are not searched for, because
   they cannot be told apart from ordinary numbers; checks 2 and 3 cover the data files.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

SCORE_COLS = [f"PC{i}" for i in range(1, 11)]
BIN_COLS = ["x_pc", "y_pc", "x_from", "x_to", "y_from", "y_to", "n"]
TOKEN_SPLIT = re.compile(r"[\s,;:\"'()\[\]{}<>=|/\\]+")


def read_text(p):
    if p.suffix == ".gz":
        with gzip.open(p, "rt", errors="replace") as f:
            return f.read()
    return p.read_text(errors="replace")


def is_num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--study", required=True)
    ap.add_argument("--share-level", required=True, choices=["individual", "aggregate"])
    ap.add_argument("--id-files", nargs="+", required=True, help="Files whose first two columns are FID IID")
    a = ap.parse_args()

    out = Path(a.out_dir)
    stem = f"{a.study}_globalPCA_v2"
    problems = []

    expected = {f"{stem}_bins.tsv", f"{stem}_coverage.tsv", f"{stem}_summary.tsv", f"{stem}.png"}
    if a.share_level == "individual":
        expected.add(f"{stem}_scores.tsv.gz")
    present = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()}
    logs = {p for p in present if re.fullmatch(r"logs/log[^/]*\.txt", p)}
    for f in sorted(expected - present):
        problems.append(f"missing expected file: {f}")
    for f in sorted(present - expected - logs):
        problems.append(f"unexpected file (move it out of {out} before upload): {f}")

    sfile = out / f"{stem}_scores.tsv.gz"
    if a.share_level == "individual" and sfile.exists():
        with gzip.open(sfile, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            if header != SCORE_COLS:
                problems.append(f"scores file columns are {header}, expected exactly {SCORE_COLS}")
            for n, line in enumerate(f, 2):
                vals = line.rstrip("\n").split("\t")
                if len(vals) != 10 or not all(is_num(v) for v in vals):
                    problems.append(f"scores file line {n} is not 10 numbers")
                    break
    bfile = out / f"{stem}_bins.tsv"
    if bfile.exists():
        lines = bfile.read_text().splitlines()
        if not lines or lines[0].split("\t") != BIN_COLS:
            problems.append(f"bins file header should be {BIN_COLS}")
        for n, line in enumerate(lines[1:], 2):
            v = line.split("\t")
            ok = (len(v) == 7 and re.fullmatch(r"PC\d+", v[0]) and re.fullmatch(r"PC\d+", v[1])
                  and all(is_num(x) for x in v[2:6]) and (v[6] == "<5" or (v[6].isdigit() and int(v[6]) >= 5)))
            if not ok:
                problems.append(f"bins file line {n} is malformed: {line}")
                break

    ids = set()
    for idf in a.id_files:
        p = Path(idf)
        if not p.exists():
            continue
        for line in p.read_text(errors="replace").splitlines():
            ids.update(line.split()[:2])
    ids.discard("0")
    searched = {i for i in ids if i and not (i.isdigit() and len(i) < 6)}
    for f in sorted(present):
        p = out / f
        if p.suffix == ".png":
            continue  # plots are drawn by project_global_pca.py and carry no sample labels
        tokens = set(TOKEN_SPLIT.split(read_text(p)))
        hits = sorted(tokens & searched)
        if hits:
            problems.append(f"{f} contains {len(hits)} sample ID(s), e.g. {hits[:3]}")

    print(f"Checked {len(present)} files in {out}; searched for {len(searched)} of {len(ids)} sample IDs "
          f"({len(ids) - len(searched)} short numeric IDs covered by the structure checks only)")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print(f"  - {p}")
        sys.exit(1)
    print("01b_pca results passed: expected files only, no sample IDs found")


if __name__ == "__main__":
    main()
