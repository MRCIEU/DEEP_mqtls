#!/usr/bin/env python3
"""DEEP global PCA v2, step 1b: convert a cohort's shared SNP file to 01b format.

Cohorts that have not run the pipeline shared SNP lists in their own formats.
This script applies the 01b SNP filters (INFO >= 0.8, MAF >= 0.01, SNP call
rate >= 0.95) and writes

    <out-root>/<cohort>/results/01/data.afreq.gz

with columns #CHROM ID REF ALT ALT_FREQS OBS_CT and canonical IDs
CHR:POS_A_B (alleles in A<C<G<T order), so that make_pca_common_snps.py
(Haotian Tang) reads it exactly like a pipeline upload.

Only autosomal biallelic A/C/G/T SNVs are kept. IDs seen more than once are
dropped entirely. GRCh37 input is lifted to GRCh38 with Hail.

Formats
  afreq    plink2 .afreq with CHR:POS_A_B IDs (e.g. Keneba, any 01b upload)
  ugr      afreq + INFO file (ID MAF R2) + post-QC SNP list (UGR)
  pelotas  plink2 .afreq with rsIDs + the matching .bim for positions (Pelotas)
  qbb      'SNPID EAF_ALL N' with IDs like 1:98929[b38]:A:G (QBB, WGS)
  actrec   Michigan FREQ csv (CHR,SNP,A1,A2,MAF,NCHROBS) + INFO csv
           (Rsq,...,CHR,BP), both sorted by CHR,BP (BMC_ACTREC)

--input is a .zip archive, a directory, or a single file. Files are located
by filename suffix, so an archive does not need unzipping (password-protected
zips must be extracted first).

Example
  python resources/global_pca/prepare_cohort_snplist.py --format qbb \
      --cohort QBB --input QBB_extracted_snp_freq.txt.zip --out-root staging/
"""

import argparse
import gzip
import io
import json
import math
import os
import re
import sys
import tempfile
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_pca_common_snps import key_id, snp_key  # noqa: E402  (canonical IDs shared with step 1c)

# Default member suffixes per format; override with --files if names differ.
DEFAULT_FILES = {
    "afreq": {"freq": ".afreq"},
    "ugr": {"freq": "_freq.afreq", "info": ".info", "keep": "_snplist.txt"},
    "pelotas": {"freq": "_freq.afreq", "bim": ".bim"},
    "qbb": {"freq": "snp_freq.txt"},
    "actrec": {"freq": "FREQ.frq_s.csv", "info": "_s.csv"},
}
ID_RE = re.compile(r"^(?:chr)?([0-9]+):([0-9]+)[_:]([ACGT]+)[_:]([ACGT]+)$")
QBB_RE = re.compile(r"^(?:chr)?([0-9XYMT]+):([0-9]+)\[b38\]:([ACGT]+):([ACGT]+)$")


# ----------------------------------------------------------------------------- input
class Source:
    """A .zip archive, a directory, or a single file; members are found by filename suffix."""

    def __init__(self, path):
        self.path = Path(path)
        self.zip = zipfile.ZipFile(self.path) if self.path.suffix == ".zip" else None
        if self.zip is None and not (self.path.is_dir() or self.path.is_file()):
            raise ValueError(f"--input not found: {self.path}")
        if self.zip:
            enc = [i.filename for i in self.zip.infolist() if i.flag_bits & 0x1]
            if enc:
                raise ValueError(f"{self.path} is password-protected; extract the needed files "
                                 f"(e.g. with 7z) into a folder and pass that folder as --input")

    def names(self):
        if self.zip:
            return [n for n in self.zip.namelist() if not n.endswith("/") and "__MACOSX" not in n]
        if self.path.is_file():
            return [str(self.path)]
        return [str(p) for p in self.path.rglob("*") if p.is_file()]

    def find(self, suffix, exclude=()):
        hits = [n for n in self.names()
                if (n.endswith(suffix) or n.endswith(suffix + ".gz"))
                and not any(n.endswith(e) or n.endswith(e + ".gz") for e in exclude)]
        if len(hits) != 1:
            raise ValueError(f"Expected one file ending '{suffix}' in {self.path}, found {hits}")
        return hits[0]

    def open(self, name):
        if self.zip:
            raw = self.zip.open(name)
            if name.endswith(".gz"):
                raw = gzip.GzipFile(fileobj=raw)
            return io.TextIOWrapper(raw, encoding="utf-8", newline=None)
        if name.endswith(".gz"):
            return gzip.open(name, "rt", newline=None)
        return open(name, "rt", newline=None)


def to_float(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return None if math.isnan(v) else v


def to_int(x):
    try:
        return int(float(x))
    except ValueError:
        return None


def limited(it, limit):
    for i, row in enumerate(it):
        if limit and i >= limit:
            return
        yield row


# --------------------------------------------------------------------------- readers
# Each reader yields (chrom, pos, a1, a2, alt_freq, obs_ct, info) with info None
# when the source has no INFO, and updates `counts` for rows it cannot parse.

def read_afreq(src, files, counts, limit, info_fail=None, keep=None):
    with src.open(files["freq"]) as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in limited(f, limit):
            r = line.rstrip("\n").split("\t")
            vid = r[ix["ID"]]
            m = ID_RE.match(vid)
            if not m:
                counts["unparsed_id"] += 1
                continue
            if m.group(1) != r[ix["#CHROM"]].removeprefix("chr"):
                counts["id_chrom_disagrees"] += 1
                continue
            if keep is not None and vid not in keep:
                counts["not_in_post_qc_snplist"] += 1
                continue
            info = None
            if info_fail is not None:
                info = 0.0 if vid in info_fail else 1.0
            yield (m.group(1), int(m.group(2)), r[ix["REF"]], r[ix["ALT"]],
                   to_float(r[ix["ALT_FREQS"]]), to_int(r[ix["OBS_CT"]]), info)


def read_ugr(src, files, counts, limit, min_info):
    # INFO and keep-list are joined by ID. Only failing IDs are held in memory.
    info_fail = set()
    with src.open(files["info"]) as f:
        header = f.readline().split()
        i_id, i_r2 = header.index("ID"), header.index("R2")
        for line in limited(f, limit):
            r = line.split()
            v = to_float(r[i_r2])
            if v is None or v < min_info:
                info_fail.add(r[i_id])
    counts["info_file_ids_below_threshold"] = len(info_fail)
    keep = None
    if "keep" in files:
        with src.open(files["keep"]) as f:
            keep = {line.strip() for line in limited(f, limit) if line.strip()}
        counts["post_qc_snplist_ids"] = len(keep)
    yield from read_afreq(src, files, counts, limit, info_fail=info_fail, keep=keep)


def read_pelotas(src, files, counts, limit):
    # The .afreq has rsIDs and no position; rows are in the same order as the .bim.
    with src.open(files["freq"]) as f, src.open(files["bim"]) as b:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for n, (line, bline) in enumerate(limited(zip(f, b), limit), 2):
            r = line.rstrip("\n").split("\t")
            bim = bline.split()
            if r[ix["ID"]] != bim[1]:
                raise ValueError(f"afreq row {n} ID {r[ix['ID']]} != bim ID {bim[1]}: files are not row aligned")
            if sorted((r[ix["REF"]], r[ix["ALT"]])) != sorted((bim[4], bim[5])):
                counts["afreq_bim_allele_mismatch"] += 1
                continue
            yield (bim[0], int(bim[3]), r[ix["REF"]], r[ix["ALT"]],
                   to_float(r[ix["ALT_FREQS"]]), to_int(r[ix["OBS_CT"]]), None)


def read_qbb(src, files, counts, limit):
    with src.open(files["freq"]) as f:
        header = f.readline().split()
        i_id, i_f, i_n = header.index("SNPID"), header.index("EAF_ALL"), header.index("N")
        for line in limited(f, limit):
            r = line.split()
            m = QBB_RE.match(r[i_id])
            if not m:
                counts["unparsed_id"] += 1
                continue
            yield (m.group(1), int(m.group(2)), m.group(3), m.group(4),
                   to_float(r[i_f]), to_int(r[i_n]), None)


def chrom_order(c):
    c = c.removeprefix("chr")
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "XY": 25, "MT": 26}.get(c, 99)


def read_actrec(src, files, counts, limit):
    """Merge-join FREQ and INFO on (CHR, BP); both must be sorted by CHR then BP.
    INFO has no alleles, so a position passes only if every INFO row there passes
    (the minimum Rsq is used)."""
    with src.open(files["freq"]) as f, src.open(files["info"]) as g:
        fh = f.readline().strip().split(",")
        gh = g.readline().strip().split(",")
        fi = {h: i for i, h in enumerate(fh)}
        rsq_col = next(h for h in gh if h.startswith("Rsq"))
        gi = {h: i for i, h in enumerate(gh)}
        info_iter = (line.strip().split(",") for line in limited(g, limit))

        def info_groups():
            # yields ((chr_order, bp), min_rsq) per position, checking sort order
            prev, cur, best = None, None, None
            for r in info_iter:
                k = (chrom_order(r[gi["CHR"]]), int(r[gi["BP"]]))
                v = to_float(r[gi[rsq_col]])
                v = -1.0 if v is None else v
                if prev is not None and k < prev:
                    raise ValueError(f"INFO file is not sorted by CHR,BP at {k}")
                prev = k
                if k != cur:
                    if cur is not None:
                        yield cur, best
                    cur, best = k, v
                else:
                    best = min(best, v)
            if cur is not None:
                yield cur, best

        groups = info_groups()
        gk, gv = next(groups, (None, None))
        prev = None
        for line in limited(f, limit):
            r = line.strip().split(",")
            parts = r[fi["SNP"]].split(":")
            if len(parts) != 4:
                # rsID rows are genotyped array SNPs; the same SNPs also appear as
                # imputed chr:pos:ref:alt rows, so nothing is lost by skipping them
                counts["skipped_rsid_rows" if r[fi["SNP"]].startswith("rs") else "unparsed_id"] += 1
                continue
            chrom, pos = r[fi["CHR"]], int(parts[1])
            k = (chrom_order(chrom), pos)
            if prev is not None and k < prev:
                raise ValueError(f"FREQ file is not sorted by CHR,BP at {k}")
            prev = k
            while gk is not None and gk < k:
                gk, gv = next(groups, (None, None))
            info = gv if gk == k else None
            if info is None:
                counts["no_info_row"] += 1
                continue
            # MAF is the frequency of A1; report A1 as ALT so ALT_FREQS = MAF
            yield (chrom, pos, r[fi["A2"]], r[fi["A1"]],
                   to_float(r[fi["MAF"]]), to_int(r[fi["NCHROBS"]]), info)


# ---------------------------------------------------------------------------- filters
def apply_filters(records, args, counts, max_obs):
    min_obs = args.min_call_rate * max_obs if max_obs else 0
    for chrom, pos, a1, a2, freq, obs, info in records:
        counts["rows_read"] += 1
        try:
            key = snp_key(chrom, pos, a1, a2)
        except ValueError:
            key = None
        if key is None:
            counts["removed_not_autosomal_biallelic_snv"] += 1
            continue
        if freq is None or obs == 0:
            counts["removed_no_observed_genotypes"] += 1
            continue
        if obs is not None and obs < min_obs:
            counts["removed_call_rate"] += 1
            continue
        if min(freq, 1 - freq) < args.min_maf:
            counts["removed_maf"] += 1
            continue
        if info is not None and info < args.min_info:
            counts["removed_info"] += 1
            continue
        counts["passed_filters"] += 1
        yield chrom.removeprefix("chr"), pos, a1, a2, freq, obs


def liftover_37_to_38(rows_path, chain, counts, threads, tmp_dir):
    """Lift (chrom, pos, a1, a2, freq, obs) rows from GRCh37 to GRCh38 with Hail."""
    import hail as hl
    hl.init(quiet=True, local=f"local[{threads}]", tmp_dir=tmp_dir, local_tmpdir=tmp_dir)
    rg37, rg38 = hl.get_reference("GRCh37"), hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(chain, rg38)
    ht = hl.import_table(rows_path, no_header=True, types={"f1": hl.tint32, "f4": hl.tfloat64, "f5": hl.tint32})
    ht = ht.annotate(lo=hl.liftover(hl.locus(ht.f0, ht.f1, "GRCh37"), "GRCh38", include_strand=True))
    ht = ht.annotate(ok=hl.is_defined(ht.lo) & (ht.lo.result.contig == "chr" + ht.f0))
    stats = ht.aggregate(hl.struct(unmapped=hl.agg.count_where(hl.is_missing(ht.lo)),
                                   other_chrom=hl.agg.count_where(hl.is_defined(ht.lo) & ~ht.ok),
                                   neg_strand=hl.agg.count_where(ht.ok & ht.lo.is_negative_strand)))
    counts["liftover_unmapped"] = stats.unmapped
    counts["liftover_other_chromosome"] = stats.other_chrom
    counts["liftover_negative_strand_alleles_complemented"] = stats.neg_strand
    ht = ht.filter(ht.ok)
    cmap = hl.literal({"A": "T", "C": "G", "G": "C", "T": "A"})
    ht = ht.select(c=ht.f0, p=ht.lo.result.position,
                   a1=hl.if_else(ht.lo.is_negative_strand, cmap[ht.f2], ht.f2),
                   a2=hl.if_else(ht.lo.is_negative_strand, cmap[ht.f3], ht.f3),
                   fq=ht.f4, ob=ht.f5)
    out = rows_path + ".b38.tsv"
    ht.key_by().select("c", "p", "a1", "a2", "fq", "ob").export(out, header=False)
    hl.stop()
    with open(out) as fh:
        for line in fh:
            c, p, a1, a2, fq, ob = line.rstrip("\n").split("\t")
            yield c, int(p), a1, a2, float(fq), int(ob)


# ------------------------------------------------------------------------------ main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--format", required=True, choices=sorted(DEFAULT_FILES))
    p.add_argument("--cohort", required=True, help="Cohort name; becomes the folder name in --out-root")
    p.add_argument("--input", required=True, help=".zip archive or directory with the cohort files")
    p.add_argument("--out-root", required=True)
    p.add_argument("--files", default=None,
                   help="Override member suffixes, e.g. 'freq=_freq.afreq,bim=.bim'")
    p.add_argument("--build", choices=["37", "38"], default="38")
    p.add_argument("--chain", default=str(Path(__file__).resolve().parents[1] /
                                          "genetics/references_grch37_to_grch38.over.chain.gz"))
    p.add_argument("--min-maf", type=float, default=0.01)
    p.add_argument("--min-info", type=float, default=0.8)
    p.add_argument("--min-call-rate", type=float, default=0.95,
                   help="Minimum OBS_CT / --n-obs-max (0 disables; rows with no observed genotypes are always dropped)")
    p.add_argument("--n-obs-max", type=int, default=None,
                   help="Allele count for full call rate (2 x n samples). If omitted, the maximum in the file is used (one extra pass).")
    p.add_argument("--limit", type=int, default=0, help="Read only the first N rows of each file (testing)")
    p.add_argument("--threads", type=int, default=4, help="Hail threads for liftover")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()

    out_dir = Path(args.out_root) / args.cohort / "results" / "01"
    out_file = out_dir / "data.afreq.gz"
    log_file = Path(args.out_root) / args.cohort / "prepare_log.json"
    if out_file.exists() and not args.overwrite:
        sys.exit(f"{out_file} exists; use --overwrite to replace it")
    out_dir.mkdir(parents=True, exist_ok=True)

    src = Source(args.input)
    files = dict(DEFAULT_FILES[args.format])
    if args.files:
        files.update(dict(kv.split("=", 1) for kv in args.files.split(",")))
    # resolve suffixes to member names (actrec: the INFO suffix also matches the FREQ file)
    resolved = {}
    for role, suffix in files.items():
        exclude = [files["freq"]] if (args.format == "actrec" and role == "info") else []
        resolved[role] = src.find(suffix, exclude=exclude)
    print(f"[{args.cohort}] files: {json.dumps(resolved)}", flush=True)

    counts = Counter()

    def records():
        if args.format == "afreq":
            return read_afreq(src, resolved, counts, args.limit)
        if args.format == "ugr":
            return read_ugr(src, resolved, counts, args.limit, args.min_info)
        if args.format == "pelotas":
            return read_pelotas(src, resolved, counts, args.limit)
        if args.format == "qbb":
            return read_qbb(src, resolved, counts, args.limit)
        return read_actrec(src, resolved, counts, args.limit)

    max_obs = args.n_obs_max
    if args.min_call_rate > 0 and max_obs is None:
        print(f"[{args.cohort}] scanning for maximum OBS_CT", flush=True)
        scan = Counter()
        max_obs = max((r[5] for r in (read_afreq(src, resolved, scan, args.limit) if args.format in ("afreq", "ugr")
                                       else read_pelotas(src, resolved, scan, args.limit) if args.format == "pelotas"
                                       else read_qbb(src, resolved, scan, args.limit) if args.format == "qbb"
                                       else read_actrec(src, resolved, scan, args.limit))
                       if r[5] is not None), default=0)
    print(f"[{args.cohort}] full-call-rate allele count: {max_obs}", flush=True)

    tmp_dir = tempfile.mkdtemp(prefix=f".{args.cohort}.", dir=out_dir)
    stage1 = os.path.join(tmp_dir, "passed.tsv")
    with open(stage1, "w") as out:
        for row in apply_filters(records(), args, counts, max_obs):
            out.write("\t".join(map(str, row)) + "\n")

    def staged_rows():
        if args.build == "37":
            yield from liftover_37_to_38(stage1, args.chain, counts, args.threads, tmp_dir)
        else:
            with open(stage1) as fh:
                for line in fh:
                    c, p_, a1, a2, fq, ob = line.rstrip("\n").split("\t")
                    yield c, int(p_), a1, a2, float(fq), int(ob)

    # drop every copy of an ID that occurs more than once
    stage2 = os.path.join(tmp_dir, "keyed.tsv")
    seen, dup = set(), set()
    with open(stage2, "w") as out:
        for c, p_, a1, a2, fq, ob in staged_rows():
            key = snp_key(c, p_, a1, a2)
            if key is None:
                counts["removed_after_liftover_not_autosomal_snv"] += 1
                continue
            if key in seen:
                dup.add(key)
            seen.add(key)
            out.write(f"{key}\t{c}\t{a1}\t{a2}\t{fq}\t{ob}\n")
    del seen
    counts["removed_duplicate_id_rows"] = 0
    n_out = 0
    with open(stage2) as fh, gzip.open(out_file, "wt", compresslevel=6) as out:
        out.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n")
        for line in fh:
            key, c, a1, a2, fq, ob = line.rstrip("\n").split("\t")
            key = int(key)
            if key in dup:
                counts["removed_duplicate_id_rows"] += 1
                continue
            out.write(f"{c}\t{key_id(key)}\t{a1}\t{a2}\t{fq}\t{ob}\n")
            n_out += 1
    counts["written"] = n_out
    for f_ in Path(tmp_dir).iterdir():
        f_.unlink()
    os.rmdir(tmp_dir)

    log = {
        "cohort": args.cohort, "format": args.format, "input": str(Path(args.input).resolve()),
        "files": resolved, "build_in": args.build, "build_out": "38",
        "filters": {"min_maf": args.min_maf, "min_info": args.min_info,
                    "min_call_rate": args.min_call_rate, "n_obs_max": max_obs},
        "limit_rows_for_testing": args.limit, "counts": dict(counts), "output": str(out_file.resolve()),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    log_file.write_text(json.dumps(log, indent=2) + "\n")
    print(f"[{args.cohort}] {json.dumps(dict(counts))}", flush=True)
    print(f"[{args.cohort}] wrote {n_out} SNPs to {out_file}", flush=True)


if __name__ == "__main__":
    main()
