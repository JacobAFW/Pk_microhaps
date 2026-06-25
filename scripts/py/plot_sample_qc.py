#!/usr/bin/env python
"""
plot_sample_qc.py — Stage 02 sanity-check plots.

Two figures:
  fws_distribution.png            — Fws histogram, threshold line, pass/fail counts.
  callable_fraction_distribution.png — callable-fraction histogram, threshold line.

Both are coloured by include_for_discovery so the cohort split is visible.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-metadata", required=True)
    p.add_argument("--fws-threshold",      type=float, required=True)
    p.add_argument("--callable-threshold", type=float, required=True)
    p.add_argument("--out-fws",      required=True)
    p.add_argument("--out-callable", required=True)
    return p.parse_args()


def load(path):
    fws, cf, include = [], [], []
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("fws") or not row.get("callable_fraction"):
                continue
            try:
                f = float(row["fws"])
                c = float(row["callable_fraction"])
            except ValueError:
                continue
            fws.append(f)
            cf.append(c)
            include.append(row["include_for_discovery"] == "TRUE")
    return np.array(fws), np.array(cf), np.array(include, dtype=bool)


def hist(values, include, threshold, xlabel, title, out_png, *, threshold_side="ge"):
    """threshold_side='ge' draws threshold and shades values >= threshold as the pass region."""
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    bins = np.linspace(min(values.min(), 0.0), max(values.max(), 1.0), 50)
    ax.hist(values[~include], bins=bins, alpha=0.6, label=f"fail ({(~include).sum()})",
            color="#c44e52")
    ax.hist(values[ include], bins=bins, alpha=0.8, label=f"pass ({include.sum()})",
            color="#4c72b0")
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1)
    ax.text(threshold, ax.get_ylim()[1] * 0.95, f" threshold = {threshold}",
            va="top", ha="left", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("samples")
    ax.set_title(title)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"==> Wrote {out_png}")


def main():
    args = parse_args()
    fws, cf, include = load(args.sample_metadata)

    if len(fws) == 0:
        sys.exit("ERROR: no usable rows in sample_metadata.tsv")

    Path(args.out_fws).parent.mkdir(parents=True, exist_ok=True)

    hist(fws, include, args.fws_threshold,
         xlabel="Fws (Pop-gen Stage 2)",
         title=f"Pk sample Fws — n={len(fws)}, pass={include.sum()}",
         out_png=args.out_fws)

    hist(cf, include, args.callable_threshold,
         xlabel="callable fraction (core+coding SNPs)",
         title=f"Pk sample callable fraction — n={len(cf)}, pass={include.sum()}",
         out_png=args.out_callable)


if __name__ == "__main__":
    main()
