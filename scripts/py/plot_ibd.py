#!/usr/bin/env python
"""
plot_ibd.py — Stage 02b sanity-check plot.

The fract_sites_IBD distribution is the health check on the whole stage. A
well-behaved cohort shows a bulk near 0 (unrelated pairs) and a small spike near
1 (clones), with little in between. Substantial mass in the middle means
residual polyclonality or a population-structure artefact from the cluster
import (decision L2) — stop and inspect rather than reading the clonal calls.

One panel per cluster plus a pooled panel, log-scaled counts so the clonal spike
stays visible against the unrelated bulk.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fract", nargs="+", required=True,
                   help="one or more <cluster>.hmm_fract.txt files")
    p.add_argument("--threshold", type=float, required=True)
    p.add_argument("--tag", default="", help="run tag, shown in the title")
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_fract(path):
    """Return fract_sites_IBD as a float array, by column name."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            col = header.index("fract_sites_IBD")
        except ValueError:
            sys.exit(f"ERROR: no fract_sites_IBD column in {path}")
        return np.array([float(ln.split("\t")[col]) for ln in fh if ln.strip()])


def main():
    args = parse_args()

    data = {}
    for path in args.fract:
        cluster = os.path.basename(path).replace(".hmm_fract.txt", "")
        data[cluster] = read_fract(path)

    pooled = np.concatenate(list(data.values())) if data else np.array([])
    panels = list(data.items()) + [("ALL CLUSTERS", pooled)]

    ncol = len(panels)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 3.8), squeeze=False)
    bins = np.linspace(0, 1, 51)

    for ax, (name, vals) in zip(axes[0], panels):
        ax.hist(vals, bins=bins, color="#4878a8", edgecolor="none")
        ax.axvline(args.threshold, color="#c44e52", linestyle="--", linewidth=1.2,
                   label=f"clone ≥ {args.threshold}")
        ax.set_yscale("log")
        n_clonal = int((vals >= args.threshold).sum())
        # Mass in the middle is the thing to look at; quantify it on the plot.
        mid = int(((vals > 0.10) & (vals < args.threshold)).sum())
        pct_mid = 100.0 * mid / len(vals) if len(vals) else 0.0
        ax.set_title(f"{name}\n{len(vals):,} pairs · {n_clonal} clonal · "
                     f"{pct_mid:.1f}% mid-range", fontsize=9)
        ax.set_xlabel("fract_sites_IBD")
        ax.set_ylabel("pairs (log)")
        ax.legend(fontsize=7, loc="upper right")

    title = "Stage 02b — IBD distribution"
    if args.tag:
        title += f"  [{args.tag}]"
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"[plot_ibd] wrote {args.out}")

    for name, vals in panels:
        mid = int(((vals > 0.10) & (vals < args.threshold)).sum())
        pct = 100.0 * mid / len(vals) if len(vals) else 0.0
        print(f"[plot_ibd] {name:<15} pairs={len(vals):<8,} "
              f"clonal={int((vals >= args.threshold).sum()):<4} "
              f"mid-range={mid} ({pct:.2f}%)")


if __name__ == "__main__":
    main()
