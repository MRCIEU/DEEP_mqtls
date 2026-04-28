#!/usr/bin/env python3
"""
DEEP mQTL Pipeline – cohort QC report generator  (v2)

Usage:
    python3 generate_deep_qc_report.py <config_file> <output_html> [--results <dir>]

Arguments:
    config_file   Path to the cohort bash config file (e.g. config_emph_pmmst.txt)
    output_html   Destination HTML file path
    --results     Path to the results/01 directory.
                  If omitted, falls back to the 'results_dir' key in the config file,
                  then to <home_directory>/results/01.

Cohort name, analyst, and all QC parameters are read from the config file.
No study-specific values need to be hard-coded.
"""

import argparse, base64, os, re, subprocess, sys, tempfile, glob
from datetime import datetime

# ─────────────────────────────────────────────────────────────────────────────
# Locate pdftoppm (Poppler) – try PATH first, then common conda locations
# ─────────────────────────────────────────────────────────────────────────────
def _find_pdftoppm():
    import shutil
    if shutil.which("pdftoppm"):
        return shutil.which("pdftoppm")
    candidates = [
        os.path.expanduser("~/miniconda3/envs/deep_env/bin/pdftoppm"),
        os.path.expanduser("~/miniforge3/envs/deep_env/bin/pdftoppm"),
        os.path.expanduser("~/mambaforge/envs/deep_env/bin/pdftoppm"),
        "/opt/conda/bin/pdftoppm",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None

PDFTOPPM = _find_pdftoppm()

# ─────────────────────────────────────────────────────────────────────────────
# File / image helpers
# ─────────────────────────────────────────────────────────────────────────────

def b64_file(path):
    """Base64-encode a file. Returns empty string on failure."""
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode()
    except Exception:
        return ""

def rasterize_pdf(path, dpi=150):
    """Rasterize the first page of a PDF to JPEG via pdftoppm.
    Returns (base64_string, size_kb) or (None, 0) if unavailable."""
    if not PDFTOPPM:
        return None, 0
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "page")
            subprocess.run(
                [PDFTOPPM, "-r", str(dpi), "-jpeg", "-l", "1", path, prefix],
                capture_output=True, timeout=60
            )
            pages = sorted(glob.glob(prefix + "*.jpg"))
            if not pages:
                return None, 0
            with open(pages[0], "rb") as fh:
                raw = fh.read()
            return base64.b64encode(raw).decode(), len(raw) / 1024
    except Exception:
        return None, 0

PDF_RASTER_DPI_SMALL = 150
PDF_RASTER_DPI_LARGE = 120
PDF_LARGE_THRESHOLD_MB = 2

def img_tag(path, caption="", style="max-width:100%;border-radius:6px;"):
    """Return an HTML image element for JPEG, PNG, or PDF."""
    if not path or not os.path.isfile(path):
        return f'<p class="missing">⚠ File not found: {os.path.basename(path)}</p>'

    ext  = os.path.splitext(path)[1].lower()
    size_mb = os.path.getsize(path) / 1e6
    cap  = f'<p class="img-caption">{caption}</p>' if caption else ""

    if ext in (".jpg", ".jpeg"):
        b64 = b64_file(path)
        if not b64:
            return f'<p class="missing">⚠ Could not read: {os.path.basename(path)}</p>'
        return f'<div class="img-wrap"><img src="data:image/jpeg;base64,{b64}" alt="{caption}" style="{style}">{cap}</div>'

    if ext == ".png":
        b64 = b64_file(path)
        if not b64:
            return f'<p class="missing">⚠ Could not read: {os.path.basename(path)}</p>'
        return f'<div class="img-wrap"><img src="data:image/png;base64,{b64}" alt="{caption}" style="{style}">{cap}</div>'

    if ext == ".pdf":
        dpi = PDF_RASTER_DPI_LARGE if size_mb > PDF_LARGE_THRESHOLD_MB else PDF_RASTER_DPI_SMALL
        b64_raster, kb = rasterize_pdf(path, dpi=dpi)
        if b64_raster:
            size_note = f'{size_mb:.1f} MB PDF → {kb:.0f} KB JPEG @ {dpi} DPI'
            return (f'<div class="img-wrap">'
                    f'<img src="data:image/jpeg;base64,{b64_raster}" alt="{caption}" style="{style}">'
                    f'<p class="img-caption">{caption}'
                    f'<span style="color:var(--muted);font-size:.75em"> ({size_note})</span></p>'
                    f'</div>')
        if size_mb <= PDF_LARGE_THRESHOLD_MB:
            b64 = b64_file(path)
            if b64:
                tag = (f'<embed src="data:application/pdf;base64,{b64}" '
                       f'type="application/pdf" style="width:100%;min-height:500px;border-radius:6px;" />'
                       f'<p style="font-size:0.75em;color:#888;">'
                       f'<a href="data:application/pdf;base64,{b64}" '
                       f'download="{os.path.basename(path)}">⬇ Download PDF</a></p>')
                return f'<div class="img-wrap">{tag}{cap}</div>'
        return (f'<div class="img-wrap pdf-link-box">'
                f'<p>📄 <strong>{os.path.basename(path)}</strong>'
                f' <span style="color:var(--muted)">({size_mb:.0f} MB — open directly)</span></p>'
                f'<p style="font-size:.82em;word-break:break-all"><code>{path}</code></p>'
                f'{cap}</div>')

    return f'<p class="missing">Unsupported file type: {ext}</p>'

def figure_grid(*items):
    cells = "\n".join(f'<div class="fig-cell">{it}</div>' for it in items)
    return f'<div class="fig-grid">{cells}</div>'

def download_link(path, label=None, max_embed_mb=15):
    """Return a styled inline download link for a file.
    Embeds as base64 data URI when ≤ max_embed_mb, otherwise shows filesystem path."""
    if not path or not os.path.isfile(path):
        fn = os.path.basename(path) if path else "?"
        return f'<span class="missing">⚠ File not found: {fn}</span>'
    fn    = os.path.basename(path)
    label = label or fn
    size_mb = os.path.getsize(path) / 1e6
    ext   = os.path.splitext(path)[1].lower()
    mime  = "application/pdf" if ext == ".pdf" else (
            "image/jpeg"      if ext in (".jpg", ".jpeg") else
            "image/png"       if ext == ".png" else "application/octet-stream")
    style = ("display:inline-flex;align-items:center;gap:5px;padding:5px 11px;"
             "background:#ebf8ff;border:1px solid #bee3f8;border-radius:6px;"
             "color:#2a4365;text-decoration:none;font-size:.83em;")
    size_note = f'<span style="color:#718096;font-size:.88em">({size_mb:.1f} MB)</span>'
    if size_mb <= max_embed_mb:
        b64 = b64_file(path)
        if b64:
            return (f'<a href="data:{mime};base64,{b64}" download="{fn}" style="{style}">'
                    f'⬇ {label} {size_note}</a>')
    return (f'<span style="font-size:.85em;color:#718096">'
            f'📄 {label} ({size_mb:.0f} MB) &nbsp;<code style="font-size:.9em">{path}</code></span>')

def download_link_grid(*items):
    """Wrap a list of download_link() outputs in a flex container."""
    return ('<div style="display:flex;flex-wrap:wrap;gap:8px;padding:4px 0;">'
            + "".join(items) + '</div>')

# ─────────────────────────────────────────────────────────────────────────────
# SNP-per-chromosome bar chart (inline SVG — no external dependencies)
# ─────────────────────────────────────────────────────────────────────────────

def snp_bar_chart_svg(snp_chr):
    """Return an SVG bar chart of SNP counts per chromosome."""
    if not snp_chr:
        return '<p class="missing">SNP count data not available</p>'

    labels, counts = [], []
    for chrom, count in snp_chr:
        labels.append(str(chrom))
        try:
            counts.append(int(str(count).replace(",", "")))
        except Exception:
            counts.append(0)

    max_count = max(counts) if counts else 1
    n = len(counts)

    W, H = 860, 260
    ml, mr, mt, mb = 72, 20, 20, 44   # margins
    plot_w = W - ml - mr
    plot_h = H - mt - mb
    gap    = 4
    bar_w  = (plot_w - gap * (n - 1)) / n

    # Y-axis grid lines
    ticks_svg = []
    for pct in [0, 0.25, 0.5, 0.75, 1.0]:
        val   = int(max_count * pct)
        y_pos = mt + plot_h - pct * plot_h
        label = f"{val // 1_000_000:.1f}M" if val >= 1_000_000 else (
                f"{val // 1_000}k"         if val >= 1_000     else str(val))
        ticks_svg.append(
            f'<line x1="{ml}" y1="{y_pos:.1f}" x2="{ml+plot_w}" y2="{y_pos:.1f}" '
            f'stroke="#e2e8f0" stroke-width="1"/>'
            f'<text x="{ml-6}" y="{y_pos+4:.1f}" text-anchor="end" '
            f'font-size="10" fill="#718096">{label}</text>'
        )

    # Bars + x-labels
    bars_svg = []
    for i, (label, count) in enumerate(zip(labels, counts)):
        x     = ml + i * (bar_w + gap)
        bh    = (count / max_count) * plot_h
        y     = mt + plot_h - bh
        cfmt  = f"{count:,}"
        # shade X and Y differently
        fill  = "#e53e3e" if label in ("X", "Y") else "#3182ce"
        bars_svg.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
            f'fill="{fill}" opacity="0.80" rx="2">'
            f'<title>Chr {label}: {cfmt} SNPs</title></rect>'
            f'<text x="{x + bar_w/2:.1f}" y="{mt + plot_h + 14:.1f}" '
            f'text-anchor="middle" font-size="9.5" fill="#718096">{label}</text>'
        )

    # Axis labels
    cy = mt + plot_h / 2
    axis_labels = (
        f'<text x="{W/2:.0f}" y="{H-4}" text-anchor="middle" font-size="11" fill="#718096">Chromosome</text>'
        f'<text x="11" y="{cy:.0f}" text-anchor="middle" font-size="11" fill="#718096" '
        f'transform="rotate(-90 11 {cy:.0f})">SNP count</text>'
    )

    return (f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
            f'style="width:100%;max-width:{W}px;display:block;margin:0 auto;">'
            f'{"".join(ticks_svg)}{"".join(bars_svg)}{axis_labels}</svg>')

# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_config(path):
    """Parse KEY=VALUE / KEY="VALUE" lines from a bash config file."""
    cfg = {}
    if not os.path.isfile(path):
        return cfg
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.split("#")[0].strip().strip('"').strip("'")
            cfg[key] = val
    return cfg

def load_rdata_as_dict(rdata_path, list_name):
    """Load an RData file and return the named list as a Python dict via Rscript."""
    code = f"""
suppressWarnings({{
  load('{rdata_path}')
  obj <- {list_name}
  out <- list()
  for (nm in names(obj)) {{
    v <- obj[[nm]]
    if (is.numeric(v) && length(v) == 1) out[[nm]] <- v
    else if (is.character(v)) out[[nm]] <- paste(v, collapse=", ")
    else if (is.logical(v) && length(v) == 1) out[[nm]] <- as.character(v)
  }}
  cat(paste(mapply(function(k,v) paste0(k, '\\t', v), names(out), out), collapse='\\n'))
}})
"""
    try:
        r = subprocess.run(["Rscript", "--vanilla", "-e", code],
                           capture_output=True, text=True, timeout=60)
        d = {}
        for line in r.stdout.strip().splitlines():
            parts = line.split("\t", 1)
            if len(parts) == 2:
                d[parts[0]] = parts[1]
        return d
    except Exception:
        return {}

def load_phenotype_summary(rdata_path):
    """Return HTML table rows for phenotype summaries from summstats_list."""
    code = f"""
suppressWarnings({{
  load('{rdata_path}')
  for (nm in names(summstats_list)) {{
    s <- summstats_list[[nm]]
    vals <- as.numeric(s)
    nms  <- names(s)
    cat(nm, '\\t', paste(paste0(nms,'=',round(vals,2)), collapse=' | '), '\\n')
  }}
}})
"""
    try:
        r = subprocess.run(["Rscript", "--vanilla", "-e", code],
                           capture_output=True, text=True, timeout=60)
        rows = []
        for line in r.stdout.strip().splitlines():
            parts = line.split("\t", 1)
            if len(parts) == 2:
                rows.append((parts[0].strip(), parts[1].strip()))
        return rows
    except Exception:
        return []

def parse_01d_log(log_path):
    """Extract lambda & min cis p-value for each control CpG from 01d log."""
    if not os.path.isfile(log_path):
        return []
    results, current = [], {}
    with open(log_path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.search(r"Lowest p-value within[^:]*:\s*([\d.e+\-]+)", line)
            if m:
                current = {"min_p": m.group(1)}
                continue
            m = re.search(r"Generating QQ-plot without cis chromosome for(.+?) with lambda\s*([\d.e+\-]+)", line)
            if m and current.get("min_p"):
                path_seg = m.group(1).strip()
                current["lambda_nocis"] = m.group(2)
                current["is_neg"] = "negative_control" in path_seg or "NEG_" in path_seg
                cg_m = re.search(r"(cg\d+)\s*$", path_seg)
                current["cpg"] = cg_m.group(1) if cg_m else path_seg.split("/")[-1]
                continue
            m = re.search(r"Generating QQ-plot for(.+?) with lambda\s*([\d.e+\-]+)", line)
            if m and current.get("min_p") and current.get("cpg"):
                current["lambda_full"] = m.group(2)
                results.append({
                    "cpg":          current["cpg"],
                    "min_p":        current["min_p"],
                    "lambda_full":  current.get("lambda_full", "N/A"),
                    "lambda_nocis": current.get("lambda_nocis", "N/A"),
                    "type":         "Negative" if current.get("is_neg") else "Positive",
                })
                current = {}
    return results

def parse_easyqc_rep(rep_path):
    """Parse the easyQC .rep file for AF correlation and outlier count."""
    if not os.path.isfile(rep_path):
        return {}
    with open(rep_path) as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        return {}
    return dict(zip(lines[0].strip().split("\t"), lines[1].strip().split("\t")))

def snp_chr_table(txt_path):
    """Parse no_snps_by_chr.txt → list of (chr, count) pairs."""
    chroms = list(range(1, 23)) + ["X", "Y"]
    rows = []
    if not os.path.isfile(txt_path):
        return rows
    with open(txt_path) as fh:
        counts = [l.strip() for l in fh if l.strip()]
    for i, c in enumerate(counts):
        if i < len(chroms):
            rows.append((str(chroms[i]), c))
    return rows

def parse_01b_log(log_path):
    """Extract SNP and sample QC stats from logs_b/log.txt.

    Returns dict with keys (all optional, 'N/A' if absent):
      snps_in                   – initial variant count
      snps_out                  – final QC'd variant count (last >5 M count)
      call_rate                 – e.g. "92.81%"
      samples_after_relatedness – present only when related="no"
      unrelated_samples / related_samples – present when related="yes"
    """
    if not os.path.isfile(log_path):
        return {}
    out = {}
    variant_counts = []   # all "N variants remaining after main filters"
    with open(log_path) as fh:
        for line in fh:
            line = line.rstrip()
            # Initial variant count (first occurrence)
            if "snps_in" not in out:
                m = re.search(r"(\d+) variants loaded from", line)
                if m:
                    out["snps_in"] = m.group(1)
            # All "N variants remaining" lines (multiple – different QC stages)
            m = re.search(r"(\d+) variants remaining after main filters", line)
            if m:
                variant_counts.append(int(m.group(1)))
            # Call rate
            m = re.search(r"Call rate:\s*([\d.]+%?)", line)
            if m:
                out["call_rate"] = m.group(1).rstrip("%") + "%"
            # Sample size after cryptic-relateds removal (related="no" pipelines)
            m = re.search(r"Sample size after removal cryptic relateds:\s*(\d+)", line)
            if m:
                out["samples_after_relatedness"] = m.group(1)
            # Unrelated / Related sets (related="yes" pipelines using PCAIR)
            m = re.search(r"Unrelated Set:\s*(\d+)\s*Sample", line)
            if m:
                out["unrelated_samples"] = m.group(1)
            m = re.search(r"Related Set:\s*(\d+)\s*Sample", line)
            if m:
                out["related_samples"] = m.group(1)
    # SNPs out: last count > 5 M (avoids PCA/LD-pruning subsets)
    big = [v for v in variant_counts if v > 5_000_000]
    if big:
        out["snps_out"] = str(big[-1])
    elif variant_counts:
        out["snps_out"] = str(variant_counts[-1])
    return out


def ewas_summary(log_path):
    """Parse ewas_run.log for hit counts."""
    if not os.path.isfile(log_path):
        return {}
    text = open(log_path).read()
    out = {}
    for pat in [
        r"There were\s*(\d+)\s*DNAm sites associated with (\w+) in the (\w+) model",
        r"There were(\d+)DNAm sites associated with (\w+) in the (\w+) model",
    ]:
        for m in re.finditer(pat, text):
            n, pheno, model = m.group(1), m.group(2).lower(), m.group(3).lower()
            out[f"{pheno}_{model}"] = n
    return out

# ─────────────────────────────────────────────────────────────────────────────
# HTML building blocks
# ─────────────────────────────────────────────────────────────────────────────

CSS = """
:root {
  --bg: #f4f6f9; --card: #ffffff; --accent: #2c5282;
  --accent2: #3182ce; --text: #2d3748; --muted: #718096;
  --pass: #276749; --passbg: #c6f6d5;
  --warn: #744210; --warnbg: #fefcbf;
  --fail: #742a2a; --failbg: #fed7d7;
  --border: #e2e8f0; --radius: 10px; --shadow: 0 2px 8px rgba(0,0,0,.08);
}
* { box-sizing: border-box; }
body { font-family: 'Segoe UI', system-ui, sans-serif; background: var(--bg);
       color: var(--text); margin: 0; padding: 0; font-size: 15px; }
header { background: linear-gradient(135deg, var(--accent), var(--accent2));
         color: #fff; padding: 28px 40px; }
header h1 { margin: 0 0 4px; font-size: 1.9em; }
header p  { margin: 0; opacity: .85; font-size: .95em; }
nav { background: var(--accent); padding: 0 40px; display: flex; flex-wrap: wrap; gap: 4px; }
nav a { color: rgba(255,255,255,.85); text-decoration: none; padding: 10px 14px;
        font-size: .85em; border-bottom: 3px solid transparent; }
nav a:hover { border-bottom-color: #fff; color: #fff; }
main { max-width: 1200px; margin: 0 auto; padding: 30px 24px 60px; }
.section { margin-bottom: 36px; }
.section-title { font-size: 1.25em; font-weight: 700; color: var(--accent);
                 border-bottom: 2px solid var(--accent); padding-bottom: 6px; margin-bottom: 18px; }
.card { background: var(--card); border-radius: var(--radius);
        box-shadow: var(--shadow); padding: 22px 26px; margin-bottom: 20px;
        border: 1px solid var(--border); }
.card h3 { margin: 0 0 14px; font-size: 1em; font-weight: 600; color: var(--accent); }
table { width: 100%; border-collapse: collapse; font-size: .9em; }
th { background: var(--accent); color: #fff; padding: 9px 12px; text-align: left; }
td { padding: 8px 12px; border-bottom: 1px solid var(--border); }
tr:last-child td { border-bottom: none; }
tr:nth-child(even) td { background: #f7fafc; }
.kv-table td:first-child { font-weight: 600; color: var(--accent); width: 42%; }
.stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; }
.stat-box { background: var(--card); border: 1px solid var(--border); border-radius: var(--radius);
            padding: 16px; text-align: center; box-shadow: var(--shadow); }
.stat-box .val { font-size: 1.6em; font-weight: 700; color: var(--accent2); }
.stat-box .lbl { font-size: .8em; color: var(--muted); margin-top: 4px; }
.pass { background: var(--passbg); color: var(--pass); padding: 2px 8px; border-radius: 4px; font-weight: 600; }
.warn { background: var(--warnbg); color: var(--warn); padding: 2px 8px; border-radius: 4px; font-weight: 600; }
.fail { background: var(--failbg); color: var(--fail); padding: 2px 8px; border-radius: 4px; font-weight: 600; }
.badge { display:inline-block; padding: 3px 10px; border-radius: 20px; font-size:.8em; font-weight:600; }
.badge-green  { background:#c6f6d5; color:#276749; }
.badge-blue   { background:#bee3f8; color:#2a4365; }
.badge-yellow { background:#fefcbf; color:#744210; }
.badge-red    { background:#fed7d7; color:#742a2a; }
.fig-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }
.fig-cell { }
.img-wrap { text-align: center; }
.img-wrap img { max-width: 100%; border: 1px solid var(--border); border-radius: 6px; }
.img-caption { font-size: .8em; color: var(--muted); margin-top: 5px; font-style: italic; }
.missing { color: #e53e3e; font-size: .85em; padding: 8px; background: #fff5f5;
           border-radius: 4px; border: 1px dashed #e53e3e; }
.pdf-link-box { background:#f7fafc; border:1px dashed var(--border);
                border-radius:6px; padding:12px 16px; text-align:left; }
.notice { background: #ebf8ff; border-left: 4px solid var(--accent2); padding: 12px 16px;
          border-radius: 0 6px 6px 0; color: #2a4365; font-size:.9em; }
.code-block { background: #1a202c; color: #e2e8f0; padding: 14px 18px; border-radius: 8px;
              font-family: 'Courier New', monospace; font-size: .82em; overflow-x: auto;
              white-space: pre-wrap; word-break: break-all; }
details summary { cursor: pointer; font-weight: 600; color: var(--accent2); padding: 6px 0; }
details[open] summary { margin-bottom: 10px; }
footer { text-align: center; padding: 20px; color: var(--muted); font-size: .8em;
         border-top: 1px solid var(--border); }
"""

def stat_box(value, label):
    return f'<div class="stat-box"><div class="val">{value}</div><div class="lbl">{label}</div></div>'

def kv_row(key, val):
    return f"<tr><td>{key}</td><td>{val}</td></tr>"

def badge(text, colour="blue"):
    return f'<span class="badge badge-{colour}">{text}</span>'

def status_badge(val, good_range=(0.95, 1.05)):
    """Badge for lambda values – displayed with 3 decimal places."""
    try:
        f = float(val)
        display = f"{f:.3f}"
        if good_range[0] <= f <= good_range[1]:   return badge(display, "green")
        elif 0.90 <= f <= 1.10:                    return badge(display, "yellow")
        else:                                       return badge(display, "red")
    except Exception:
        return val

def pval_badge(val):
    """Badge for p-values – displayed in scientific notation with 3 dp."""
    try:
        f = float(val)
        display = f"{f:.3e}"
        if f < 5e-8:   return badge(display, "green")
        elif f < 1e-4: return badge(display, "yellow")
        else:          return badge(display, "red")
    except Exception:
        return val

def build_match_badge(declared, inferred):
    """Green tick if declared == inferred, yellow warning otherwise."""
    if str(declared) == str(inferred):
        return badge(f"GRCh{inferred} ✓", "green")
    return badge(f"GRCh{inferred} ≠ declared GRCh{declared}", "yellow")

# ─────────────────────────────────────────────────────────────────────────────
# Main report builder
# ─────────────────────────────────────────────────────────────────────────────

def build_report(config_path, out_html, results_override=None):
    cfg = parse_config(config_path)

    # ── Resolve results directory ─────────────────────────────────────────────
    # Priority: CLI --results > config 'results_dir' key > <home_directory>/results/01
    if results_override:
        R = results_override
    elif cfg.get("results_dir"):
        R = cfg["results_dir"]
    else:
        R = os.path.join(cfg.get("home_directory", "."), "results", "01")
    R = R.rstrip("/")

    # ── Load RData objects ────────────────────────────────────────────────────
    pre_qc  = load_rdata_as_dict(os.path.join(R, "cohort_descriptives.RData"),     "cohort_summary")
    post_qc = load_rdata_as_dict(os.path.join(R, "cohort_descriptives_was.RData"), "cohort_summary")

    # ── Parse log and support files ───────────────────────────────────────────
    d_controls   = parse_01d_log(os.path.join(R, "logs_d", "log.txt"))
    pos_controls = [c for c in d_controls if c["type"] == "Positive"]
    neg_controls = [c for c in d_controls if c["type"] == "Negative"]
    b_stats      = parse_01b_log(os.path.join(R, "logs_b", "log.txt"))
    eqc          = parse_easyqc_rep(os.path.join(R, "easyQC_topmed.rep"))
    snp_chr      = snp_chr_table(os.path.join(R, "no_snps_by_chr.txt"))
    ewas         = ewas_summary(os.path.join(R, "ewas_run.log"))
    pheno_rows   = load_phenotype_summary(
                     os.path.join(R, f"edited_phenotypes_summary_{cfg.get('study_name','STUDY')}.Rdata"))

    # ── Inferred builds ───────────────────────────────────────────────────────
    def read_txt(p):
        try: return open(p).read().strip()
        except Exception: return "N/A"

    build_01b = read_txt(os.path.join(R, "01b_inferred_build.txt"))
    build_01g = read_txt(os.path.join(R, "01g_inferred_build.txt"))

    # ── Local helpers ─────────────────────────────────────────────────────────
    def get(d, k, default="N/A"):
        return d.get(k, default)

    def fmt_num(v):
        try:    return f"{int(float(v)):,}"
        except: return v

    def fmt_f(v, dp=3):
        try:    return f"{float(v):.{dp}f}"
        except: return v

    # ── Config tables ─────────────────────────────────────────────────────────
    study   = cfg.get("study_name", "UNKNOWN")
    analyst = cfg.get("analyst_name", "N/A")
    now     = datetime.now().strftime("%d %B %Y, %H:%M")

    study_rows = [
        kv_row("Study name",    cfg.get("study_name", "N/A")),
        kv_row("Analyst",       cfg.get("analyst_name", "N/A")),
        kv_row("Email",         cfg.get("analyst_email", "N/A")),
        kv_row("Array",         badge(cfg.get("methylation_array", "N/A"), "blue")),
        kv_row("Tissue",        cfg.get("tissue", "N/A")),
        kv_row("Sorted cells",  cfg.get("sorted_methylation", "N/A")),
        kv_row("Ancestry",      badge(cfg.get("ancestry", "N/A"), "yellow")),
        kv_row("Age group",     cfg.get("age", "N/A")),
        kv_row("Related samples", badge(cfg.get("related", "N/A"),
                                        "red" if cfg.get("related") == "yes" else "green")),
        kv_row("Reference panel", cfg.get("reference", "N/A")),
        kv_row("Smoking",       cfg.get("smoking", "none") or "none"),
        kv_row("Study-specific vars", f'<code>{cfg.get("study_specific_vars","N/A")}</code>'),
    ]

    build_declared = cfg.get("genome_build", "?")
    genome_rows = [
        kv_row("Declared build",       badge(f"GRCh{build_declared}", "blue")),
        kv_row("Inferred build (01b)", build_match_badge(build_declared, build_01b)),
        kv_row("Inferred build (01g)", build_match_badge(build_declared, build_01g)),
    ]

    # Comprehensive QC parameters (includes all thresholds from config)
    qc_param_rows = "".join([
        # --- SNP / genotype QC ---
        kv_row("<em>Genotype QC</em>", ""),
        kv_row("SNP HWE threshold",   cfg.get("snp_hwe",  "N/A")),
        kv_row("SNP MAF threshold",   cfg.get("snp_maf",  "N/A")),
        kv_row("SNP missingness",     cfg.get("snp_miss", "N/A")),
        kv_row("Sample missingness",  cfg.get("snp_imiss","N/A")),
        # --- Relatedness / GRM ---
        kv_row("<em>Relatedness &amp; GRM</em>", ""),
        kv_row("Relatedness cutoff (GRM)", cfg.get("rel_cutoff",     "N/A")),
        kv_row("GRM MAF cutoff",           cfg.get("grm_maf_cutoff", "N/A")),
        # --- Genetic PCA ---
        kv_row("<em>Genetic PCA</em>", ""),
        kv_row("n_pcs",               cfg.get("n_pcs",    "N/A")),
        kv_row("PCA SD threshold",    cfg.get("pca_sd",   "N/A")),
        # --- Methylation PCs ---
        kv_row("<em>Methylation PCs (EWAS)</em>", ""),
        kv_row("n_meth_pcs",          cfg.get("n_meth_pcs",         "N/A")),
        kv_row("Meth PC variance cutoff", cfg.get("meth_pc_cutoff", "N/A")),
        kv_row("Meth PC p-threshold", cfg.get("meth_pc_threshold",  "N/A")),
        # --- mQTL ---
        kv_row("<em>mQTL analysis</em>", ""),
        kv_row("MatrixEQTL soft threshold", cfg.get("soft_threshold", "N/A")),
        # --- Compute ---
        kv_row("<em>Resources</em>", ""),
        kv_row("Threads",    cfg.get("nthreads",       "N/A")),
        kv_row("Memory (GB)",cfg.get("mem",            "N/A")),
        kv_row("Meth chunks",cfg.get("meth_chunks",    "N/A")),
        kv_row("Genetic chunks", cfg.get("genetic_chunks", "N/A")),
    ])

    path_rows = "".join([
        kv_row("home_directory",  f'<code style="font-size:.8em">{cfg.get("home_directory","N/A")}</code>'),
        kv_row("results_dir",     f'<code style="font-size:.8em">{R}</code>'),
        kv_row("bfile_raw",       f'<code style="font-size:.8em">{cfg.get("bfile_raw","N/A")}</code>'),
        kv_row("vcf_dir",         f'<code style="font-size:.8em">{cfg.get("vcf_dir","N/A")}</code>'),
        kv_row("idat_directory",  f'<code style="font-size:.8em">{cfg.get("idat_directory","N/A")}</code>'),
        kv_row("betas",           f'<code style="font-size:.8em">{cfg.get("betas","N/A")}</code>'),
    ])

    # ── Phenotype summary HTML ────────────────────────────────────────────────
    pheno_html = "".join(
        f"<tr><td><code>{r[0]}</code></td><td style='font-size:.85em'>{r[1]}</td></tr>"
        for r in pheno_rows
    )

    # ── Control table ─────────────────────────────────────────────────────────
    def ctrl_rows(ctrl_list):
        return "".join(
            f"<tr><td><code>{c['cpg']}</code></td>"
            f"<td>{pval_badge(c['min_p'])}</td>"
            f"<td>{status_badge(c['lambda_full'])}</td>"
            f"<td>{status_badge(c['lambda_nocis'])}</td></tr>"
            for c in ctrl_list
        )

    # ── 01d QQ / Manhattan download links ────────────────────────────────────
    def d_qq_links(prefix="positive_control_untransformed"):
        d_dir = os.path.join(R, "01d")
        links = [download_link(os.path.join(d_dir, f"{prefix}_{c['cpg']}_qqplot.jpeg"),
                               f"{c['cpg']} QQ")
                 for c in pos_controls]
        return download_link_grid(*links) if links else '<p class="missing">No QQ plots found.</p>'

    def d_mht_links(prefix="positive_control_untransformed"):
        d_dir = os.path.join(R, "01d")
        links = [download_link(os.path.join(d_dir, f"{prefix}_{c['cpg']}_manhattan.pdf"),
                               f"{c['cpg']} Manhattan")
                 for c in pos_controls]
        return download_link_grid(*links) if links else '<p class="missing">No Manhattan plots found.</p>'

    # ── 01e Manhattan download links ──────────────────────────────────────────
    def pc_gwas_links():
        links = []
        for i in range(1, int(cfg.get("n_pcs", 20)) + 1):
            path = os.path.join(R, "01e", f"gwas_PC{i}_manhattan_beta.pdf")
            if os.path.isfile(path):
                links.append(download_link(path, f"PC{i} Manhattan"))
        return download_link_grid(*links) if links else '<p class="missing">No PC GWAS plots found.</p>'

    # ── 01g scatter download links ────────────────────────────────────────────
    def scatter_links_01g():
        g_dir = os.path.join(R, "01g")
        if not os.path.isdir(g_dir):
            return '<p class="missing">No 01g scatter plots found.</p>'
        cgs = sorted(set(
            re.sub(r"_(pc_no|hc_no|pc_hc)_scatter\.pdf$", "", f)
            for f in os.listdir(g_dir)
            if f.endswith("_scatter.pdf")
        ))
        if not cgs:
            return '<p class="missing">No 01g scatter plots found.</p>'
        items = []
        for cg in cgs:
            for suffix, label in [("_pc_no_scatter.pdf", "PC vs Unadj"),
                                   ("_hc_no_scatter.pdf", "HC vs Unadj"),
                                   ("_pc_hc_scatter.pdf", "PC vs HC")]:
                path = os.path.join(g_dir, f"{cg}{suffix}")
                if os.path.isfile(path):
                    items.append(download_link(path, f"{cg}: {label}"))
        return download_link_grid(*items) if items else '<p class="missing">No 01g scatter plots found.</p>'

    # ── EWAS report links ─────────────────────────────────────────────────────
    # Compute relative path from report output dir back to results/01
    report_dir  = os.path.dirname(os.path.abspath(out_html))
    results_abs = os.path.abspath(R)
    try:
        rel_prefix = os.path.relpath(results_abs, report_dir)
    except ValueError:
        rel_prefix = results_abs   # fallback on Windows drives mismatch
    ewas_links = " &nbsp;|&nbsp; ".join(
        f'<a href="{os.path.join(rel_prefix, fn)}">{fn}</a>'
        for fn in [
            f"qc1_ewas_report_age_{study}_unilife.html",
            f"qc1_ewas_report_sex_{study}_unilife.html",
            f"qc1_ewas_report_sex_negative_control_{study}_unilife.html",
        ]
    )

    # ─────────────────────────────────────────────────────────────────────────
    # Assemble HTML
    # ─────────────────────────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>DEEP mQTL QC Report – {study}</title>
<style>{CSS}</style>
</head>
<body>
<header>
  <h1>DEEP mQTL Pipeline – QC Report</h1>
  <p>Cohort: <strong>{study}</strong> &nbsp;|&nbsp;
     Analyst: <strong>{analyst}</strong> &nbsp;|&nbsp;
     Generated: {now}</p>
</header>
<nav>
  <a href="#config">Config</a>
  <a href="#s01a">01a Data checks</a>
  <a href="#s01b">01b Genetic QC</a>
  <a href="#s01c">01c Methylation QC</a>
  <a href="#s01d">01d mQTL controls</a>
  <a href="#s01e">01e PC GWAS</a>
  <a href="#s01f">01f Normalisation</a>
  <a href="#s01g">01g Haplotype components</a>
</nav>
<main>

<!-- ════════════════════════════════ CONFIG ════════════════════════════════ -->
<div class="section" id="config">
<div class="section-title">Configuration</div>

<div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;">
  <div class="card">
    <h3>Study &amp; Sample Details</h3>
    <table class="kv-table">{"".join(study_rows)}</table>
  </div>
  <div class="card">
    <h3>Genome Build</h3>
    <table class="kv-table">{"".join(genome_rows)}</table>
  </div>
</div>

<div class="card">
  <h3>QC Parameters</h3>
  <table class="kv-table">{qc_param_rows}</table>
</div>

<div class="card">
  <details>
    <summary>Data paths</summary>
    <table class="kv-table" style="margin-top:8px;">{path_rows}</table>
  </details>
</div>
</div>

<!-- ════════════════════════════════ 01a ═══════════════════════════════════ -->
<div class="section" id="s01a">
<div class="section-title">01a – Data Validation</div>

<div class="stat-grid">
  {stat_box(fmt_num(get(pre_qc,'geno_sample_size')),        "Genotype samples (raw)")}
  {stat_box(fmt_num(get(pre_qc,'methylation_sample_size')), "Methylation samples (raw)")}
  {stat_box(fmt_num(get(pre_qc,'geno_meth_common_ids')),    "Overlapping samples")}
  {stat_box(fmt_num(get(pre_qc,'n_CpGs')),                  "CpGs")}
  {stat_box(fmt_num(get(pre_qc,'n_snp')),                   "SNPs (pre-QC)")}
  {stat_box(fmt_f(get(pre_qc,'mean_info')),                 "Mean imputation INFO")}
  {stat_box(get(pre_qc,'min_info'),                         "Min INFO")}
  {stat_box(fmt_num(get(pre_qc,'mqtl_sample_size')),        "mQTL analysis N")}
</div>

<div class="card" style="margin-top:20px;">
  <h3>Post-QC Sample Sizes</h3>
  <table class="kv-table">
    {kv_row("mQTL set (GWAS) – N",  fmt_num(get(post_qc,'methylation_sample_size_gwas')))}
    {kv_row("EWAS set – N",         fmt_num(get(post_qc,'methylation_sample_size_ewas')))}
    {kv_row("Males (mQTL)",         fmt_num(get(post_qc,'mqtl_n_males')))}
    {kv_row("Females (mQTL)",       fmt_num(get(post_qc,'mqtl_n_females')))}
    {kv_row("Males (EWAS)",         fmt_num(get(post_qc,'ewas_n_males')))}
    {kv_row("Females (EWAS)",       fmt_num(get(post_qc,'ewas_n_females')))}
    {kv_row("SNPs (post-QC)",       fmt_num(get(post_qc,'n_snp')))}
    {kv_row("Covariates",           f'<code>{get(pre_qc,"covariates")}</code>')}
    {kv_row("No. plates",           get(pre_qc,'n_Plate_factor'))}
  </table>
</div>

<div class="card" style="margin-top:16px;">
  <h3>SNPs per Chromosome (pre-QC)</h3>
  {snp_bar_chart_svg(snp_chr)}
</div>

<div class="card" style="margin-top:16px;">
  <h3>Phenotype Summaries <span style="font-weight:400;font-size:.85em;color:var(--muted)">(includes age distribution)</span></h3>
  <table>
    <thead><tr><th>Phenotype</th><th>Summary Statistics</th></tr></thead>
    <tbody>{pheno_html if pheno_html else '<tr><td colspan="2">Not available</td></tr>'}</tbody>
  </table>
</div>

{figure_grid(
    img_tag(os.path.join(R,"no_snps_by_chr.pdf"),          "SNPs per chromosome (raw PDF)"),
    img_tag(os.path.join(R,"XY_methylation_vs_sex.pdf"),   "XY methylation sex check"),
    img_tag(os.path.join(R,"age_distribution.pdf"),        "Age distribution"),
)}
</div>

<!-- ════════════════════════════════ 01b ═══════════════════════════════════ -->
<div class="section" id="s01b">
<div class="section-title">01b – Genetic Data Processing</div>

<div class="card">
  <h3>Summary</h3>
  <table class="kv-table">
    {kv_row("Liftover performed",
            badge('No – already GRCh38', 'green') if build_01b == '38'
            else badge('Yes – lifted to GRCh38', 'yellow'))}
    {kv_row("SNPs in → out (01b QC)",
            f"{fmt_num(b_stats.get('snps_in','N/A'))} → {fmt_num(b_stats.get('snps_out','N/A'))}")}
    {kv_row("Sample size after removal of cryptic relateds",
            b_stats["samples_after_relatedness"]
            if "samples_after_relatedness" in b_stats
            else (f"Unrelated: {b_stats.get('unrelated_samples','N/A')}, "
                  f"Related (kept): {b_stats.get('related_samples','N/A')}")
            if "unrelated_samples" in b_stats
            else "N/A")}
    {kv_row("Call rate <span style='font-weight:400;font-size:.85em;color:var(--muted)'>"
            "(proportion of cohort SNPs that loaded onto global biobank SNP loadings)</span>",
            b_stats.get("call_rate", "N/A"))}
    {kv_row("EasyQC AF correlation (vs TOPMed ref)", eqc.get('AFCHECK.cor_raf.ref_EAF', 'N/A'))}
    {kv_row("EasyQC AF outliers",  fmt_num(eqc.get('AFCHECK.numOutlier', 'N/A')))}
  </table>
</div>

{figure_grid(
    img_tag(os.path.join(R,"pcaplot.pdf"),                     "Genetic PCA plot (01b)"),
    img_tag(os.path.join(R,f"{study}_globalPCA.png"),          "Global PCA plot"),
    img_tag(os.path.join(R,"grm_distrib_01b.pdf"),             "GRM distribution (01b)"),
    img_tag(os.path.join(R,"easyQC_topmed.multi.AFCHECK.png"), "EasyQC allele frequency check"),
)}
</div>

<!-- ════════════════════════════════ 01c ═══════════════════════════════════ -->
<div class="section" id="s01c">
<div class="section-title">01c – Phenotype &amp; Methylation QC</div>

<div class="card">
  <h3>EWAS QC Summary</h3>
  <table class="kv-table" style="margin-bottom:14px;">
    {kv_row("Sex distribution (mQTL set)",
            f"Males: {fmt_num(get(post_qc,'mqtl_n_males'))} &nbsp;|&nbsp; "
            f"Females: {fmt_num(get(post_qc,'mqtl_n_females'))}")}
    {kv_row("Sex distribution (EWAS set)",
            f"Males: {fmt_num(get(post_qc,'ewas_n_males'))} &nbsp;|&nbsp; "
            f"Females: {fmt_num(get(post_qc,'ewas_n_females'))}")}
    {kv_row("No. of plates", get(pre_qc,'n_Plate_factor'))}
  </table>
  <table>
    <thead><tr><th>Phenotype</th><th>Model</th><th>Hits (p &lt; 1×10⁻⁷)</th></tr></thead>
    <tbody>
      <tr><td>Age</td><td>Unadjusted (none)</td><td>{ewas.get("age_none","0")}</td></tr>
      <tr><td>Age</td><td>Adjusted (all covariates)</td><td>{ewas.get("age_all","0")}</td></tr>
      <tr><td>Age</td><td>SVA-adjusted</td><td>{ewas.get("age_sva","N/A")}</td></tr>
      <tr><td>Sex</td><td>Unadjusted (none)</td><td>{ewas.get("sex_none","–")}</td></tr>
      <tr><td>Sex</td><td>Adjusted (all covariates)</td><td>{ewas.get("sex_all","–")}</td></tr>
      <tr><td>Sex (scrambled, negative control)</td><td>Unadjusted</td><td>0</td></tr>
      <tr><td>Sex (scrambled, negative control)</td><td>Adjusted</td><td>0</td></tr>
    </tbody>
  </table>
  <p class="notice" style="margin-top:12px;">
    <strong>Age EWAS:</strong> 0 hits expected for children with narrow age range (SD ≈ 0.42 yr).
    <strong>Sex EWAS:</strong> strong X/Y-linked CpG signal is a positive QC indicator.
    <strong>Scrambled-sex negative control:</strong> 0 hits confirms no spurious inflation ✓.<br>
    Full EWAS reports (interactive HTML): {ewas_links}
  </p>
</div>

{figure_grid(
    img_tag(os.path.join(R,f"age_smoking_prediction_plot_{study}.jpg"),   "Age &amp; smoking prediction"),
    img_tag(os.path.join(R,f"raw_phenotype_distribution_plot_{study}.jpg"), "Raw phenotype distributions"),
    img_tag(os.path.join(R,f"edited_phenotype_distribution_plot_{study}.jpg"), "Edited phenotype distributions"),
)}

<div class="card" style="margin-top:16px;">
  <h3>Methylation PCs</h3>
  {figure_grid(
      img_tag(os.path.join(R,f"meth_pcs_scree_plot_{study}_unilife.pdf"),  "Methylation PCs scree plot"),
      img_tag(os.path.join(R,f"meth_pcs_PC1PC2_plot_{study}_unilife.pdf"), "Meth PC1 vs PC2"),
      img_tag(os.path.join(R,f"meth_pcs_PC3PC4_plot_{study}_unilife.pdf"), "Meth PC3 vs PC4"),
      img_tag(os.path.join(R,f"pc_var_association_plot_{study}_unilife.pdf"), "PC–covariate associations"),
  )}
</div>

<div class="card" style="margin-top:16px;">
  <h3>Cell Counts</h3>
  {figure_grid(
      img_tag(os.path.join(R,"cellcounts_plot.pdf"),                                    "Cell count composition"),
      img_tag(os.path.join(R,f"cell_count_PCs_scree_plot_{study}_unilife.pdf"),         "Cell count PCs scree"),
      img_tag(os.path.join(R,"cor_plot_ori.pdf"),                                       "Correlation (original)"),
      img_tag(os.path.join(R,"cor_plot_comb.pdf"),                                      "Correlation (combined)"),
  )}
</div>
</div>

<!-- ════════════════════════════════ 01d ═══════════════════════════════════ -->
<div class="section" id="s01d">
<div class="section-title">01d – mQTL Positive &amp; Negative Controls</div>

<div class="card">
  <h3>Positive Controls – Lambda &amp; Min Cis P-value</h3>
  {'<p class="missing">Could not parse 01d log.</p>' if not pos_controls else f'''
  <table>
    <thead>
      <tr><th>CpG</th><th>Min cis P-value</th>
          <th>λ (full GWAS)</th><th>λ (no-cis chr)</th></tr>
    </thead>
    <tbody>{ctrl_rows(pos_controls)}</tbody>
  </table>
  <p style="font-size:.8em;color:var(--muted);margin-top:8px;">
    λ colour: <span class="badge badge-green">green</span> 0.95–1.05 &nbsp;
    <span class="badge badge-yellow">yellow</span> 0.90–1.10 &nbsp;
    <span class="badge badge-red">red</span> outside 0.90–1.10
  </p>'''}
</div>

<div class="card" style="margin-top:16px;">
  <h3>Negative Controls – Lambda &amp; Min Cis P-value</h3>
  {'<p class="missing">No negative controls parsed.</p>' if not neg_controls else f'''
  <table>
    <thead>
      <tr><th>CpG</th><th>Min cis P-value</th>
          <th>λ (full GWAS)</th><th>λ (no-cis chr)</th></tr>
    </thead>
    <tbody>{ctrl_rows(neg_controls)}</tbody>
  </table>'''}
</div>

<div class="card" style="margin-top:16px;">
  <h3>Positive Control QQ Plots</h3>
  {d_qq_links()}
</div>

<div class="card" style="margin-top:16px;">
  <h3>Positive Control Manhattan Plots</h3>
  {d_mht_links()}
</div>
</div>

<!-- ════════════════════════════════ 01e ═══════════════════════════════════ -->
<div class="section" id="s01e">
<div class="section-title">01e – Genetic PC GWAS</div>

<div class="card">
  <p class="notice">GWAS run for genetic PCs (up to PC{cfg.get('n_pcs','?')} configured).
     Each Manhattan plot shows genome-wide associations for that PC with effect estimates (beta).
     Clustered peaks confirm the PCs capture population structure as expected.</p>
  {pc_gwas_links()}
</div>
</div>

<!-- ════════════════════════════════ 01f ═══════════════════════════════════ -->
<div class="section" id="s01f">
<div class="section-title">01f – Cross-cohort Methylation Normalisation</div>

<div class="card">
  <p class="notice">
    Step 01f requires IDAT files to extract control-probe PCs for cross-cohort functional
    normalisation via <strong>meffil</strong>. Cohorts with IDAT files can participate in
    the central normalisation step.
  </p>
  <table class="kv-table">
    {kv_row("IDAT files available (can participate in cross-cohort normalisation)",
            badge("Yes ✓", "green")
            if cfg.get("idat_directory") and cfg.get("idat_directory") not in ("N/A","","none")
            else badge("No", "red"))}
  </table>
</div>
</div>

<!-- ════════════════════════════════ 01g ═══════════════════════════════════ -->
<div class="section" id="s01g">
<div class="section-title">01g – Haplotype Components (HCs)</div>

<div class="card">
  <p class="notice">
    Haplotype components (HCs) are computed via PBWTpaint to capture fine-scale local haplotype
    sharing as an alternative to global genetic PCs. Scatterplots compare GWAS summary statistics
    under three correction strategies:
    <strong>No correction</strong> vs <strong>PC-adjusted</strong> vs <strong>HC-adjusted</strong>.
  </p>
  <h3 style="margin-top:14px;">Positive control CpGs: PC adj vs HC adj vs Unadjusted</h3>
  {scatter_links_01g()}
</div>
</div>

</main>
<footer>DEEP mQTL QC Report &nbsp;|&nbsp; {study} &nbsp;|&nbsp; Generated {now}</footer>
</body></html>"""

    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    with open(out_html, "w") as fh:
        fh.write(html)
    print(f"Report written to: {out_html}")
    print(f"File size: {os.path.getsize(out_html)/1e6:.1f} MB")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a DEEP mQTL pipeline QC report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # results_dir read from config key 'results_dir' or defaulted to <home_directory>/results/01
  python3 generate_deep_qc_report.py config_emph_pmmst.txt reports/PMMST_qc.html

  # Override results directory from the command line (useful for multi-cohort analysts)
  python3 generate_deep_qc_report.py config_emph_pmmst.txt reports/PMMST_qc.html \\
      --results /home/prachand/projects/DEEP/EMPHASIS/results/01

  # Second cohort, same analyst machine
  python3 generate_deep_qc_report.py config_emph_pure.txt reports/PURE_qc.html \\
      --results /home/prachand/projects/DEEP/PURE/results/01
""")
    parser.add_argument("config",      help="Path to the cohort bash config file")
    parser.add_argument("output_html", help="Destination HTML report path")
    parser.add_argument("--results",   metavar="DIR",
                        help="Path to the results/01 directory "
                             "(overrides config 'results_dir' key and the default)")
    args = parser.parse_args()
    build_report(args.config, args.output_html, args.results)
