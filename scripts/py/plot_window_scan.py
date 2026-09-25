#!/usr/bin/env python
"""
plot_window_scan.py — Stage 04 figures.

Four figures, each written as a vector PDF and a 300-dpi PNG preview
(Nature Genetics standard, per MISSION's output standard):

  window_het_distribution     He histogram over all unique windows, 0.5 line
                              marked.                     (NB1 c23 / NB2 c23)
  window_het_vs_entropy       He-vs-entropy scatter + He boxplot by variant
                              count. This is THE direct visual check against
                              the paper: the boxplot must rise monotonically
                              and flatten past ~10 variants.          (NB2 c33)
  window_manhattan_het        Per-window He along all 14 chromosomes, with
                              candidate windows highlighted.     (NB2 plot_man)
  window_het_raw_vs_nomiss    het_raw vs het_nomiss coloured by n_variants —
                              quantifies deviation D1.                  (ours)

Reads `windows_all.tsv.gz` (all unique non-empty windows) and takes the
heterozygome thresholds as arguments so the figures annotate the same cut the
scan applied.

Usage:
    python scripts/py/plot_window_scan.py \
        --windows outputs/scan/v1_2026-09-25/windows_all.tsv.gz \
        --min-snps 3 --min-het 0.50 \
        --outdir reports/figures/04_window_scan
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DPI = 300
FIGURES = ("window_het_distribution", "window_het_vs_entropy",
           "window_manhattan_het", "window_het_raw_vs_nomiss")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--windows", required=True, help="windows_all.tsv.gz from window_scan.py")
    p.add_argument("--min-snps", type=int, required=True)
    p.add_argument("--min-het", type=float, required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def save(fig, outdir, stem):
    """Vector PDF plus 300-dpi PNG preview, same stem."""
    pdf = Path(outdir) / f"{stem}.pdf"
    png = Path(outdir) / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {pdf} + {png}")


# ---------------------------------------------------------------------------
def fig_het_distribution(df, min_het, outdir):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(df["het_raw"].dropna(), bins=60, color="#4C72B0", edgecolor="white",
            linewidth=0.3, label="het_raw")
    ax.hist(df["het_nomiss"].dropna(), bins=60, histtype="step", color="#C44E52",
            linewidth=1.2, label="het_nomiss")
    ax.axvline(min_het, color="black", linestyle="--", linewidth=1)
    ax.annotate(f"scan floor He = {min_het:g}", xy=(min_het, ax.get_ylim()[1] * 0.92),
                xytext=(6, 0), textcoords="offset points", fontsize=8)
    ax.set_xlabel("window heterozygosity (He = 1 - $\\Sigma f^2$)")
    ax.set_ylabel("windows")
    ax.set_title(f"Stage 04 — window He, all {len(df):,} unique windows")
    ax.legend(frameon=False, fontsize=8)
    save(fig, outdir, "window_het_distribution")


def fig_het_vs_entropy(df, min_het, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    sc = ax.scatter(df["het_raw"], df["entropy"], c=df["n_variants"], s=3,
                    cmap="viridis", alpha=0.45, rasterized=True)
    ax.axvline(min_het, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("het_raw")
    ax.set_ylabel("entropy ($-\\Sigma f \\ln f$)")
    ax.set_title("He vs entropy")
    fig.colorbar(sc, ax=ax, label="variants in window")

    # He by variant count — the shape check against Siegel: monotone rise,
    # diminishing returns past ~10 variants.
    ax = axes[1]
    counts = sorted(df["n_variants"].unique())
    shown = [c for c in counts if c <= 15]
    data = [df.loc[df["n_variants"] == c, "het_raw"].dropna().values for c in shown]
    tail = df.loc[df["n_variants"] > 15, "het_raw"].dropna().values
    labels = [str(c) for c in shown]
    if tail.size:
        data.append(tail)
        labels.append(">15")
    ax.boxplot(data, labels=labels, showfliers=False,
               medianprops=dict(color="#C44E52"))
    ax.axhline(min_het, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("variants in window")
    ax.set_ylabel("het_raw")
    ax.set_title("He by variant count")
    fig.tight_layout()
    save(fig, outdir, "window_het_vs_entropy")


def fig_manhattan(df, min_snps, min_het, outdir):
    chroms = sorted(df["chrom"].unique())
    fig, ax = plt.subplots(figsize=(13, 4.5))
    offset, ticks, tick_labels = 0, [], []
    cand_mask = (df["n_variants"] >= min_snps) & (df["het_raw"] >= min_het)
    for i, c in enumerate(chroms):
        sub = df[df["chrom"] == c]
        x = sub["midpoint"].values + offset
        base = "#9AA7B8" if i % 2 == 0 else "#C9D2DC"
        ax.scatter(x, sub["het_raw"], s=1.2, color=base, alpha=0.6, rasterized=True)
        csub = sub[cand_mask.loc[sub.index]]
        ax.scatter(csub["midpoint"].values + offset, csub["het_raw"], s=1.8,
                   color="#C44E52" if i % 2 == 0 else "#DD8452", alpha=0.8,
                   rasterized=True)
        span = sub["window_end"].max()
        ticks.append(offset + span / 2)
        tick_labels.append(c.replace("ordered_PKNH_", "").replace("_v2", ""))
        offset += span
    ax.axhline(min_het, color="black", linestyle="--", linewidth=0.8)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlim(0, offset)
    ax.set_xlabel("chromosome")
    ax.set_ylabel("het_raw")
    ax.set_title(f"Stage 04 — per-window He across 14 chromosomes "
                 f"(coloured: {int(cand_mask.sum()):,} candidate windows, "
                 f"$\\geq${min_snps} SNPs & He $\\geq$ {min_het:g})")
    save(fig, outdir, "window_manhattan_het")


def fig_raw_vs_nomiss(df, outdir):
    sub = df.dropna(subset=["het_nomiss"])
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(sub["het_nomiss"], sub["het_raw"], c=sub["n_variants"], s=3,
                    cmap="viridis", alpha=0.45, rasterized=True)
    lim = [0, max(1.0, float(sub[["het_raw", "het_nomiss"]].max().max()))]
    ax.plot(lim, lim, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("het_nomiss (missing-bearing haplotypes dropped)")
    ax.set_ylabel("het_raw (Siegel-faithful)")
    med = float((sub["het_raw"] - sub["het_nomiss"]).median())
    ax.set_title(f"Deviation D1 — median $\\Delta$He = {med:.4f}")
    fig.colorbar(sc, ax=ax, label="variants in window")
    save(fig, outdir, "window_het_raw_vs_nomiss")


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.windows, sep="\t")
    print(f"==> {len(df):,} unique windows from {args.windows}")

    fig_het_distribution(df, args.min_het, outdir)
    fig_het_vs_entropy(df, args.min_het, outdir)
    fig_manhattan(df, args.min_snps, args.min_het, outdir)
    fig_raw_vs_nomiss(df, outdir)

    missing = [s for s in FIGURES
               if not ((outdir / f"{s}.pdf").exists() and (outdir / f"{s}.png").exists())]
    if missing:
        raise SystemExit(f"ERROR: figures not written: {missing}")
    print(f"==> 4 figures x (PDF + {DPI} dpi PNG) in {outdir}")


if __name__ == "__main__":
    main()
