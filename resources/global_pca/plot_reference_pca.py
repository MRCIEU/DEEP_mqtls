#!/usr/bin/env python3
"""DEEP global PCA v2, step 3: reference PCA figures (central use, not run by cohorts).

Writes
  <out>_reference_check.png  PC1/PC2, PC1/PC3, PC2/PC4; one column per region highlighted
  <out>_reference_3D.png     PC1-PC3 from two fixed viewing angles
  <out>_reference_3D.html    interactive PC1-PC3 (needs plotly; skipped if not installed)

Example
  python resources/global_pca/plot_reference_pca.py \
      --scores resources/global_pca/reference/deep_globalpca_v2_scores_flat.tsv.gz \
      --eigenvalues resources/global_pca/reference/deep_globalpca_v2_eigenvalues.txt \
      --out plots/deep_globalpca_v2
"""

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REGIONS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]
NAMES = {"AFR": "Africa", "AMR": "Americas", "CSA": "Central & South Asia", "EAS": "East Asia",
         "EUR": "Europe", "MID": "Middle East", "OCE": "Oceania"}
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
INK, INK2, MUTED, GRID, CTX, HI = "#1f1f1d", "#5f5e58", "#8a8983", "#e6e5e0", "#d6d5cf", "#2a78d6"
PAIRS = [(1, 2), (1, 3), (2, 4)]


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": GRID,
                         "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5})


def check_grid(d, pct, n_snps, out):
    fig, axes = plt.subplots(3, 8, figsize=(19, 7.6), constrained_layout=True)
    for r, (a, b) in enumerate(PAIRS):
        x, y = f"PC{a}", f"PC{b}"
        rx, ry = d[x].max() - d[x].min(), d[y].max() - d[y].min()
        lim = ((d[x].min() - rx * .06, d[x].max() + rx * .06), (d[y].min() - ry * .06, d[y].max() + ry * .06))
        for c in range(8):
            ax = axes[r, c]
            if c == 0:
                ax.scatter(d[x], d[y], s=3, c="#9c9b95", linewidths=0, alpha=0.6, rasterized=True)
                placed = []
                for reg in sorted(REGIONS, key=lambda g: -(d.genetic_region == g).sum()):
                    s = d[d.genetic_region == reg]
                    if s.empty:
                        continue
                    px, py = s[x].median(), s[y].median()
                    for _ in range(12):  # nudge labels apart
                        if all(abs((px - qx) / rx) > 0.13 or abs((py - qy) / ry) > 0.075 for qx, qy in placed):
                            break
                        py += ry * 0.075 if py >= np.median([q[1] for q in placed]) else -ry * 0.075
                    placed.append((px, py))
                    ax.text(px, py, reg, fontsize=7.5, color=INK, ha="center", va="center", fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
                if r == 0:
                    ax.set_title(f"All reference samples\n(n = {len(d):,})", fontsize=8.5, color=INK)
            else:
                reg = REGIONS[c - 1]
                m = d.genetic_region == reg
                ax.scatter(d.loc[~m, x], d.loc[~m, y], s=3, c=CTX, linewidths=0, rasterized=True)
                ax.scatter(d.loc[m, x], d.loc[m, y], s=5, c=HI, linewidths=0, alpha=0.85, rasterized=True)
                if r == 0:
                    ax.set_title(f"{NAMES[reg]}\n(n = {m.sum():,})", fontsize=8.5, color=INK)
            ax.set_xlim(lim[0]); ax.set_ylim(lim[1])
            ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(length=2)
            if c == 0:
                ax.set_ylabel(f"PC{b} ({pct[b]:.1f}%)", fontsize=8)
            else:
                ax.set_yticklabels([])
            ax.set_xlabel(f"PC{a} ({pct[a]:.1f}%)", fontsize=7.5)
    fig.suptitle(f"DEEP global PCA v2 reference: HGDP + 1000 Genomes, {len(d):,} unrelated samples, {n_snps:,} SNPs. "
                 "% = share of variance among the top 20 PCs", fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.savefig(f"{out}_reference_check.png", dpi=150, facecolor="white")
    plt.close(fig)


def three_d_png(d, n_snps, out):
    cent = d.groupby("genetic_region")[["PC1", "PC2", "PC3"]].median().reindex(REGIONS).dropna()
    fig = plt.figure(figsize=(13.5, 6.6))
    views = [(25, 150, "View A", {"EUR": (-0.02, -0.03, -0.03), "MID": (0.02, 0.03, 0.03)}),
             (40, -120, "View B", {"AMR": (-0.04, 0.0, 0.20), "OCE": (0.06, 0.0, 0.06), "CSA": (0.07, 0.0, 0.0),
                                   "MID": (-0.04, 0.03, 0.0), "EUR": (0, -0.05, -0.04)})]
    for i, (e, az, title, off) in enumerate(views):
        a = fig.add_axes([0.0 + i * 0.52, 0.1, 0.47, 0.83], projection="3d")
        for reg, col in zip(REGIONS, COLOURS):
            s = d[d.genetic_region == reg]
            a.scatter(s.PC1, s.PC2, s.PC3, s=3.5, c=col, alpha=0.85, linewidths=0, depthshade=False,
                      label=f"{NAMES[reg]} ({len(s):,})")
        for reg in cent.index:
            dx, dy, dz = off.get(reg, (0, 0, 0.045))
            a.text(cent.loc[reg, "PC1"] + dx, cent.loc[reg, "PC2"] + dy, cent.loc[reg, "PC3"] + dz, reg, fontsize=8.5,
                   fontweight="bold", color=INK, ha="center",
                   bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
        a.view_init(elev=e, azim=az)
        a.set_box_aspect(None, zoom=1.08)
        for lab, fn in (("PC1", a.set_xlabel), ("PC2", a.set_ylabel), ("PC3", a.set_zlabel)):
            fn(lab, color=INK2, labelpad=1)
        for axis in (a.xaxis, a.yaxis, a.zaxis):
            axis.set_pane_color((1, 1, 1, 1))
            axis._axinfo["grid"]["color"] = GRID
            axis._axinfo["grid"]["linewidth"] = 0.5
        a.tick_params(labelsize=6.5, colors=MUTED, pad=0)
        a.set_title(title, color=INK, fontsize=9, y=0.98)
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=7, frameon=False, markerscale=3, fontsize=8.5,
               columnspacing=2.2, handletextpad=0.4)
    fig.text(0.01, 0.97, f"DEEP global PCA v2 reference, PC1 to PC3 (HGDP + 1000 Genomes, {len(d):,} samples, "
             f"{n_snps:,} SNPs)", fontsize=10.5, color=INK)
    fig.savefig(f"{out}_reference_3D.png", dpi=160, facecolor="white")
    plt.close(fig)


def three_d_html(d, n_snps, out):
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not installed: skipping the interactive 3D plot")
        return
    fig = go.Figure()
    for reg, col in zip(REGIONS, COLOURS):
        s = d[d.genetic_region == reg]
        fig.add_trace(go.Scatter3d(
            x=s.PC1, y=s.PC2, z=s.PC3, mode="markers", name=f"{NAMES[reg]} ({len(s):,})",
            marker=dict(size=2.6, color=col, opacity=0.85, line=dict(width=0)),
            customdata=np.stack([s.s, s.population], axis=1),
            hovertemplate="<b>%{customdata[1]}</b> (" + reg + ")<br>%{customdata[0]}<br>"
                          "PC1 %{x:.3f}<br>PC2 %{y:.3f}<br>PC3 %{z:.3f}<extra></extra>"))
    cent = d.groupby("genetic_region")[["PC1", "PC2", "PC3"]].median().reindex(REGIONS).dropna()
    fig.add_trace(go.Scatter3d(x=cent.PC1, y=cent.PC2, z=cent.PC3 + 0.03, mode="text", text=list(cent.index),
                               textfont=dict(size=13, color=INK), showlegend=False, hoverinfo="skip"))
    ax = dict(backgroundcolor="white", gridcolor=GRID, zerolinecolor="#c9c8c2", showspikes=False,
              tickfont=dict(size=10, color=MUTED), title_font=dict(size=13, color=INK2))
    r, e, az = 2.1, np.radians(40), np.radians(-120)
    eye = dict(x=r * np.cos(e) * np.cos(az), y=r * np.cos(e) * np.sin(az), z=r * np.sin(e))
    fig.update_layout(
        title=dict(text=f"DEEP global PCA v2 reference, PC1 to PC3<br><sup>HGDP + 1000 Genomes, {len(d):,} samples, "
                        f"{n_snps:,} SNPs. Drag to rotate, scroll to zoom, click a region in the legend to hide or "
                        "show it.</sup>", x=0.01, font=dict(size=16, color=INK)),
        scene=dict(xaxis=dict(title="PC1", **ax), yaxis=dict(title="PC2", **ax), zaxis=dict(title="PC3", **ax),
                   camera=dict(eye=eye), aspectmode="cube"),
        legend=dict(title=dict(text="Reference region (n)"), itemsizing="constant", font=dict(size=12, color=INK),
                    x=0.01, y=0.85),
        paper_bgcolor="white", margin=dict(l=0, r=0, t=70, b=0), hoverlabel=dict(bgcolor="white", font_size=12))
    fig.write_html(f"{out}_reference_3D.html", include_plotlyjs=True, full_html=True,
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 3}})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scores", required=True, help="<prefix>_scores_flat.tsv.gz from make_reference_pca.py")
    p.add_argument("--eigenvalues", required=True, help="<prefix>_eigenvalues.txt")
    p.add_argument("--snplist", default=None, help="<prefix>.snplist (for the SNP count in titles)")
    p.add_argument("--out", required=True, help="Output prefix")
    args = p.parse_args()
    d = pd.read_csv(args.scores, sep="\t")
    ev = pd.read_csv(args.eigenvalues, sep="\t")
    pct = {int(r.PC[2:]): 100 * r.prop_of_top_k for r in ev.itertuples()}
    n_snps = sum(1 for _ in open(args.snplist)) if args.snplist else 0
    style()
    check_grid(d, pct, n_snps, args.out)
    three_d_png(d, n_snps, args.out)
    three_d_html(d, n_snps, args.out)
    print(f"Wrote {args.out}_reference_check.png, _reference_3D.png and _reference_3D.html")


if __name__ == "__main__":
    main()
