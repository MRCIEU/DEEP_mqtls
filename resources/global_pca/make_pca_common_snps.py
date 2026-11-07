#!/usr/bin/env python3
"""Extract raw GRCh38 reference SNP IDs and intersect with cohort 01b afreqs.

Only autosomal biallelic A/C/G/T SNVs are eligible. No MAF, missingness,
FILTER, sample, LD, or PCA filtering is applied. Allele order is ignored;
complementary strands are NOT merged. Cohort inputs must already be GRCh38.
Outputs are gzip, one canonical CHR:POS_A_B ID per line, without a header.
The JSON manifest is published last, after validation. Existing outputs are
never overwritten. Use a new output prefix for a new input snapshot.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone

BASES = 'ACGT'


def chromosome_code(value):
    value = value.removeprefix('chr')
    return {'X': '23', 'Y': '24', 'XY': '25', 'MT': '26', 'M': '26'}.get(value, value)


def snp_key(chrom, pos, a, b):
    chrom = str(chrom).removeprefix('chr')
    if not chrom.isdigit() or not 1 <= int(chrom) <= 22:
        return None
    pos = int(pos)
    if not 0 < pos < 2**32:
        raise ValueError('Invalid position: %s' % pos)
    if len(a) != 1 or len(b) != 1 or a not in BASES or b not in BASES or a == b:
        return None
    a, b = sorted((BASES.index(a), BASES.index(b)))
    return ((int(chrom) * 2**32 + pos) << 4) | (a * 4 + b)


def key_id(key):
    locus, alleles = divmod(key, 16)
    chrom, pos = divmod(locus, 2**32)
    return '%d:%d_%s_%s' % (chrom, pos, BASES[alleles // 4], BASES[alleles % 4])


def id_key(identifier):
    locus, a, b = identifier.split('_')
    chrom, pos = locus.split(':')
    return snp_key(chrom, pos, a, b)


def fingerprint(path):
    path = Path(path).resolve()
    st = path.stat()
    return {'path': str(path), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns}


def cohort_inputs(roots):
    chosen = {}
    skipped = []
    for root in map(Path, roots):
        if not root.is_dir():
            raise ValueError('Missing cohort root: %s' % root)
        for folder in sorted(root.iterdir()):
            if not folder.is_dir():
                continue
            path = folder / 'results/01/data.afreq.gz'
            if folder.name in chosen:
                skipped.append({'cohort': folder.name, 'directory': str(folder), 'reason': 'preferred_input_already_selected'})
            elif path.is_file():
                chosen[folder.name] = path
            else:
                skipped.append({'cohort': folder.name, 'directory': str(folder), 'reason': 'no_afreq'})
    if not chosen:
        raise ValueError('No cohort SNP lists found')
    return chosen, skipped


def read_cohort(path):
    keys = set()
    counts = Counter()
    with gzip.open(path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        indices = [header.index(x) for x in ('#CHROM', 'ID', 'REF', 'ALT')]
        for lineno, row in enumerate(reader, 2):
            counts['input_rows'] += 1
            try:
                chrom, identifier, ref, alt = [row[i] for i in indices]
                locus, a, b = identifier.split('_')
                id_chrom, pos = locus.split(':')
                if chromosome_code(chrom) != chromosome_code(id_chrom):
                    raise ValueError('ID chromosome disagrees with #CHROM')
                if sorted((a, b)) != sorted((ref, alt)):
                    raise ValueError('ID alleles disagree with REF/ALT')
                key = snp_key(chrom, pos, ref, alt)
                if key is None:
                    counts['excluded_non_autosomal_biallelic_snv'] += 1
                    continue
                if key in keys:
                    raise ValueError('Duplicate canonical SNP ID: %s' % identifier)
                keys.add(key)
            except (ValueError, IndexError) as exc:
                raise ValueError('%s:%d: %s' % (path, lineno, exc)) from exc
    if not keys:
        raise ValueError('No eligible SNPs: %s' % path)
    counts['eligible_unique_snps'] = len(keys)
    return keys, dict(counts)


def reference_keys(lines, chromosome, counts):
    previous_pos = -1
    group = set()
    for line in lines:
        chrom, pos, a, b = line.rstrip('\n').split('\t')
        counts['input_variant_records'] += 1
        if chrom.removeprefix('chr') != str(chromosome):
            raise ValueError('Unexpected reference chromosome: %s' % chrom)
        pos = int(pos)
        if pos < previous_pos:
            raise ValueError('Reference is not coordinate sorted')
        if pos != previous_pos:
            yield from sorted(group)
            group = set()
            previous_pos = pos
        key = snp_key(chrom, pos, a, b)
        if key is None:
            counts['excluded_non_biallelic_snv'] += 1
            continue
        if key in group:
            raise ValueError('Duplicate reference SNP: %s' % key_id(key))
        group.add(key)
        counts['raw_snps'] += 1
    yield from sorted(group)


def extract_reference(paths, bcftools, output, logdir):
    stats = []
    with gzip.open(output, 'wt', compresslevel=1) as out:
        for chromosome, path in enumerate(paths, 1):
            counts = Counter()
            with (logdir / ('bcftools.chr%d.log' % chromosome)).open('w') as err:
                command = [bcftools, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', str(path)]
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True)
                try:
                    for key in reference_keys(proc.stdout, chromosome, counts):
                        out.write(key_id(key) + '\n')
                    proc.stdout.close()
                    if proc.wait() != 0:
                        raise RuntimeError('bcftools query failed: %s; see %s' % (path, err.name))
                finally:
                    if proc.poll() is None:
                        proc.terminate()
                        proc.wait()
            if not counts['raw_snps']:
                raise ValueError('No SNPs on chromosome %d' % chromosome)
            stats.append({'chromosome': chromosome, **dict(counts)})
            print('Reference chr%d: %s' % (chromosome, dict(counts)), flush=True)
    return stats


def validate_lists(raw, common, expected_raw, expected_common):
    # Single merge scan checks gzip integrity, canonical IDs, strict order,
    # uniqueness, counts, and that every common SNP appears in raw.
    with gzip.open(raw, 'rt') as r, gzip.open(common, 'rt') as c:
        target = next(c, None)
        count = matched = 0
        previous = -1
        for line in r:
            key = id_key(line.strip())
            if key is None or key <= previous or line.strip() != key_id(key):
                raise ValueError('Invalid raw SNP list order/ID')
            previous = key
            count += 1
            if line == target:
                matched += 1
                target = next(c, None)
        if target is not None or count != expected_raw or matched != expected_common:
            raise ValueError('Output counts/subset validation failed')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bcf-dir', type=Path, required=True)
    p.add_argument('--pipeline-qc-dir', required=True)
    p.add_argument('--pipeline-debug-dir', required=True)
    p.add_argument('--bcftools', required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--prefix', default='hgdp_1kg.GRCh38')
    p.add_argument('--raw-snplist', type=Path, help='Reuse and preserve an existing raw reference SNP list.')
    p.add_argument('--resume-stage', type=Path,
                   help='Reuse completed raw extraction in a failed staging directory; input fingerprints must still match.')
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / (args.prefix + '.lock')).open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('Another run with prefix %s holds %s; wait for it or use a new prefix' % (args.prefix, lock.name))
        except OSError as exc:
            # DEEP global PCA v2 change: NFS mounts (e.g. /mnt/sod) may not support flock
            print('WARNING: could not lock %s (%s); continuing without a lock' % (lock.name, exc), flush=True)
        files = {k: args.output_dir / (args.prefix + suffix) for k, suffix in {
            'raw': '.raw.biallelic_snps.snplist.gz',
            'common': '.all_cohorts_common.snplist.gz',
            'cohort_common': '.cohorts_only_common.snplist.gz',
            'cohorts': '.cohort_inputs.tsv', 'manifest': '.manifest.json'}.items()}
        if any(path.exists() for key, path in files.items()
               if not (key == 'raw' and args.raw_snplist)):
            raise FileExistsError('Outputs for prefix %s already exist in %s; use a new prefix' % (args.prefix, args.output_dir))
        cohorts, skipped = cohort_inputs([args.pipeline_qc_dir, args.pipeline_debug_dir])
        paths = [args.bcf_dir / ('hgdp1kgp_chr%d.filtered.SNV_INDEL.phased.shapeit5.bcf' % c) for c in range(1, 23)]
        snapshots = [fingerprint(x) for x in paths + list(cohorts.values())]
        tool_version = subprocess.check_output([args.bcftools, '--version'], text=True).splitlines()[0]
        raw_snapshot = fingerprint(args.raw_snplist) if args.raw_snplist else None
        print('Selected cohorts (%d): %s' % (len(cohorts), ', '.join(cohorts)), flush=True)
        # Staging stays alongside final files, never in /tmp. Keep failed work
        # and stderr logs for diagnosis; successful work removes staged files.
        stage = args.resume_stage or Path(tempfile.mkdtemp(prefix='.' + args.prefix + '.', dir=args.output_dir))
        if stage.resolve().parent != args.output_dir.resolve():
            raise ValueError('Staging directory must be directly inside output directory')
        staged = {k: stage / path.name for k, path in files.items()}
        checkpoint_path = stage / 'raw_checkpoint.json'
        if args.raw_snplist:
            if args.resume_stage:
                raise ValueError('Use raw-snplist or resume-stage, not both')
            staged['raw'] = args.raw_snplist.resolve()
            files['raw'] = staged['raw']
            raw_stats = []
        elif args.resume_stage:
            checkpoint = json.loads(checkpoint_path.read_text())
            if (checkpoint['reference_inputs'] != snapshots[:22]
                    or checkpoint['raw_file'] != fingerprint(staged['raw'])
                    or checkpoint['bcftools'] != tool_version
                    or checkpoint['raw_schema'] != 1):
                raise ValueError('Raw extraction checkpoint does not match current inputs/tool')
            raw_stats = checkpoint['reference_counts']
            if [r['chromosome'] for r in raw_stats] != list(range(1, 23)):
                raise ValueError('Incomplete raw extraction checkpoint')
            print('Resuming validated-input raw extraction checkpoint', flush=True)
        else:
            raw_stats = extract_reference(paths, args.bcftools, staged['raw'], stage)
            if snapshots[:22] != [fingerprint(x) for x in paths]:
                raise RuntimeError('Reference inputs changed during extraction')
            checkpoint_path.write_text(json.dumps({
                'raw_schema': 1, 'reference_inputs': snapshots[:22],
                'reference_counts': raw_stats, 'raw_file': fingerprint(staged['raw']),
                'bcftools': tool_version}, indent=2) + '\n')
        common = None
        cohort_stats = []
        # Smaller inputs first reduce the in-memory intersection quickly.
        for name, path in sorted(cohorts.items(), key=lambda item: (item[1].stat().st_size, item[0])):
            keys, counts = read_cohort(path)
            if common is None:
                common = keys
            else:
                common.intersection_update(keys)
                del keys
            cohort_stats.append({'cohort': name, **fingerprint(path), **counts,
                                 'remaining_cohort_intersection': len(common)})
            print('Cohort %s: %s; intersection=%d' % (name, counts, len(common)), flush=True)
            if not common:
                raise ValueError('Empty intersection after %s' % name)
        pre_ref_count = len(common)
        with gzip.open(staged['cohort_common'], 'wt', compresslevel=1) as out:
            for key in sorted(common):
                out.write(key_id(key) + '\n')
        with gzip.open(staged['cohort_common'], 'rt') as src:
            if [id_key(line.strip()) for line in src] != sorted(common):
                raise ValueError('Cohort-only intersection validation failed')
        os.replace(staged['cohort_common'], files['cohort_common'])
        print('Saved cohort-only intersection: %d SNPs' % pre_ref_count, flush=True)
        retained = 0
        scanned_raw = 0
        with gzip.open(staged['raw'], 'rt') as src, gzip.open(staged['common'], 'wt', compresslevel=1) as dst:
            for line in src:
                scanned_raw += 1
                key = id_key(line.strip())
                if key in common:
                    dst.write(line)
                    common.remove(key)
                    retained += 1
        if not retained:
            raise ValueError('No cohort-common SNPs match reference')
        raw_count = scanned_raw if args.raw_snplist else sum(x['raw_snps'] for x in raw_stats)
        print('Validating raw=%d and common=%d SNPs' % (raw_count, retained), flush=True)
        validate_lists(staged['raw'], staged['common'], raw_count, retained)
        if snapshots != [fingerprint(x) for x in paths + list(cohorts.values())]:
            raise RuntimeError('Inputs changed during the build')
        if raw_snapshot and raw_snapshot != fingerprint(args.raw_snplist):
            raise RuntimeError('Raw SNP list changed during the build')
        with staged['cohorts'].open('w') as out:
            # DEEP global PCA v2 change: union of keys, so cohorts with different count fields can be written
            fields = list(dict.fromkeys(k for row in cohort_stats for k in row))
            writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', restval=0)
            writer.writeheader()
            writer.writerows(cohort_stats)
        manifest = {
            'status': 'validated', 'completed_utc': datetime.now(timezone.utc).isoformat(),
            'genome_build': 'GRCh38 (input contract: official phased BCF and normalized DEEP 01b afreq)',
            'definition': 'Autosomal biallelic A/C/G/T SNVs; unordered alleles; no strand complementation; no FILTER/MAF/genotype/sample/LD filtering',
            'raw_snps': raw_count, 'cohort_common_before_reference': pre_ref_count,
            'common_snps': retained, 'cohort_common_absent_from_reference': len(common),
            'bcftools': tool_version, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'reference_inputs': snapshots[:22], 'reference_counts': raw_stats,
            'reused_raw_input': raw_snapshot,
            'cohorts': cohort_stats, 'skipped_directories': skipped,
            'outputs': {k: str(v) for k, v in files.items()},
        }
        staged['manifest'].write_text(json.dumps(manifest, indent=2) + '\n')
        for key in ('raw', 'common', 'cohorts', 'manifest'):
            if key == 'raw' and args.raw_snplist:
                continue
            os.replace(staged[key], files[key])
        if checkpoint_path.exists():
            checkpoint_path.unlink()
        # Preserve nonempty bcftools warnings; discard empty stderr files.
        for log in stage.glob('bcftools.*.log'):
            if log.stat().st_size == 0:
                log.unlink()
        if not list(stage.iterdir()):
            stage.rmdir()
        print('Validated outputs: %s' % json.dumps(manifest['outputs']), flush=True)


if __name__ == '__main__':
    main()
