#!/usr/bin/env python
"""
plot_snp_qc.py — Stage 03 sanity-check plots.

Two figures:
  snp_qc_maf_density.png         — MAF density on the kept SNPs (left-truncated at maf_min).
  snp_qc_missingness_density.png — F_MISSING density (right-truncated at missingness_max).

Reads MAF + F_MISSING from the bcftools-tagged VCF via `bcftools query`.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf",              required=True)
    p.add_argument("--maf-min",          type=float, required=True)
    p.add_argument("--missingness-max",  type=float, required=True)
    p.add_argument("--out-maf",          required=True)
    p.add_argument("--out-missingness",  required=True)
    return p.parse_args()


def bcftools_query(vcf, fmt):
    """Stream the requested INFO fields from the VCF."""
    cmd = ["bcftools", "query", "-f", fmt, vcf]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        sys.exit("ERROR: bcftools not on PATH. Activate the vvg-box env first.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR: bcftools query failed:\n{e.stderr}")
    return res.stdout.splitlines()


def parse_floats(lines):
    out = []
    for ln in lines:
        ln = ln.strip()
        if not ln or ln == ".":
            continue
        try:
            out.append(float(ln))
        except ValueError:
            continue
    return np.array(out)


def density(values, threshold, *, side, xlabel, title, out_png):
    """side='lower' shades pass region as values >= threshold (MAF);
       side='upper' shades pass region as values <= threshold (missingness)."""
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    if len(values) == 0:
        ax.text(0.5, 0.5, "no SNPs", ha="center", va="center", transform=ax.transAxes)
    else:
        bins = np.linspace(0, max(values.max(), 1.0), 60)
        ax.hist(values, bins=bins, alpha=0.85, color="#4c72b0",
                label=f"kept SNPs (n={len(values)})")
        ax.axvline(threshold, color="black", linestyle="--", linewidth=1)
        ax.text(threshold, ax.get_ylim()[1] * 0.95, f" threshold = {threshold}",
                va="top", ha="left", fontsize=9)
        # Sanity: with side='lower' all values should be >= threshold, etc.
        if side == "lower":
            n_violate = int((values < threshold).sum())
        else:
            n_violate = int((values > threshold).sum())
        if n_violate:
            ax.text(0.98, 0.02, f"WARNING: {n_violate} SNPs violate the threshold",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8, color="#c44e52")
        ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("SNPs")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"==> Wrote {out_png}")


def main():
    args = parse_args()

    # MAF = min(AF, 1-AF). We requested AF in build_discovery_snp_vcf.
    af_lines = bcftools_query(args.vcf, "%INFO/AF\n")
    af = parse_floats(af_lines)
    maf = np.minimum(af, 1.0 - af) if len(af) else af

    miss_lines = bcftools_query(args.vcf, "%INFO/F_MISSING\n")
    miss = parse_floats(miss_lines)

    Path(args.out_maf).parent.mkdir(parents=True, exist_ok=True)

    density(maf, args.maf_min, side="lower",
            xlabel="MAF",
            title=f"Pk discovery cohort — MAF density (n={len(maf)} SNPs)",
            out_png=args.out_maf)

    density(miss, args.missingness_max, side="upper",
            xlabel="F_MISSING",
            title=f"Pk discovery cohort — per-SNP missingness (n={len(miss)} SNPs)",
            out_png=args.out_missingness)


if __name__ == "__main__":
    main()
