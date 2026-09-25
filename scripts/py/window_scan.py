#!/usr/bin/env python
"""
window_scan.py — Stage 04: the Pk heterozygome (Siegel NB1 phase 3).

Slides a fixed-width window along the discovery SNP set one chromosome at a
time, scores every window that contains at least one SNP, and writes the
candidate-window table ("heterozygome") that Stage 05 selects the panel from.

Ported from `external/clones/vivax-mhaps/notebooks/1_evaluating_marker_candidates.ipynb`
cells 36-41. The notebook's markdown says windows with ">1 variant" are stored; its
code keeps `n_variants != 0`. This follows the CODE.

No dask (deviation D3): 593 samples x 51,251 variants is ~60 MB as int8.

---------------------------------------------------------------------------
LOCKED DECISIONS (handoffs/04_window_scan.md; confirmed at phase-0 review)
---------------------------------------------------------------------------
L1  Missing-call policy. Both are emitted on every row of both outputs:
      het_raw     Siegel-faithful — haplotypes containing a missing call are
                  counted as distinct alleles (only flagged, never dropped).
      het_nomiss  haplotypes containing any missing call dropped before
                  frequencies are computed.
    The heterozygome is filtered on het_raw ONLY, so the scan stays comparable
    to the paper. delta_het = het_raw - het_nomiss travels with every row.

L2  Het-call policy. het_raw / entropy / n_alleles / n_alleles_with_missing /
    n_alleles_with_het are computed on the UN-HAPLOIDIFIED diploid sample
    columns, exactly as NB1 cell 36 does. haploidify_samples() with Siegel's
    seed 250523 is applied ONLY to the het_nomiss branch. Haploidifying before
    het_raw would break comparability with the paper and make delta_het
    meaningless. The seed is written into the summary.

L3  Window placement vs core boundaries. frac_in_core (fraction of the window
    span overlapping pk_core_coding.bed) is REPORTED, never filtered on.
    Filtering is Stage 05's call.

L4  Window de-duplication. Windows carrying an identical set of variant
    positions are collapsed, keeping the first occurrence — NB1's
    `np.unique(positions, return_index=True)` block. This materially changes
    the window count, so both pre- and post-dedup counts go in the summary.
---------------------------------------------------------------------------

Usage:
    python scripts/py/window_scan.py \
        --vcf outputs/qc/snps.discovery.vcf.gz \
        --bed data/reference/pk_core_coding.bed \
        --samples outputs/qc/discovery_samples.declonal.txt \
        --snps-summary outputs/qc/snps.summary.tsv \
        --window-size 200 --step 50 --min-snps 3 --min-het 0.50 --seed 250523 \
        --out-windows outputs/scan/v1_2026-09-25/windows_all.tsv.gz \
        --out-heterozygome outputs/scan/v1_2026-09-25/heterozygome.tsv \
        --out-summary outputs/scan/v1_2026-09-25/window_scan_summary.tsv

Guards (all fail loud, none silent):
  * no output path may already exist — this stage never overwrites
  * VCF sample count must equal the Stage 02b cohort list
  * variants read must equal outputs/qc/snps.summary.tsv, overall AND per chromosome
  * the L4 de-duplication function carries its own self-test

Exit codes:
    0  scan completed
    1  a guard failed (message on stderr says which)
"""

import argparse
import sys
from pathlib import Path

import allel
import numpy as np
import pandas as pd

# Column order is fixed by the handoff's output contract. Do not reorder.
COLUMNS = [
    "chrom", "window_start", "window_end", "midpoint",
    "n_variants", "variant_positions",
    "n_alleles", "het_raw", "entropy",
    "n_alleles_with_missing", "n_alleles_with_het",
    "n_alleles_nomiss", "het_nomiss", "entropy_nomiss",
    "delta_het", "n_samples_used", "frac_in_core",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Stage 03 discovery SNP VCF")
    p.add_argument("--bed", required=True, help="pk_core_coding.bed (for frac_in_core, L3)")
    p.add_argument("--samples", required=True,
                   help="Stage 02b cohort list; asserted against the VCF header")
    p.add_argument("--snps-summary", required=True,
                   help="Stage 03 snps.summary.tsv; asserted against variants read")
    p.add_argument("--window-size", type=int, required=True)
    p.add_argument("--step", type=int, required=True)
    p.add_argument("--min-snps", type=int, required=True)
    p.add_argument("--min-het", type=float, required=True)
    p.add_argument("--seed", type=int, required=True, help="haploidify seed (L2)")
    p.add_argument("--missing-policy", required=True,
                   help="L1. Only 'both' is implemented — het_raw and het_nomiss are "
                        "always both emitted and the heterozygome is always filtered "
                        "on het_raw. Asserted, not branched on: a value other than "
                        "'both' is a config error, not a silent fallback.")
    p.add_argument("--out-windows", required=True)
    p.add_argument("--out-heterozygome", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------
def refuse_to_overwrite(paths):
    """Protocol standing rule: nothing existing gets overwritten."""
    clashes = [p for p in paths if Path(p).exists()]
    if clashes:
        sys.exit("ERROR: these outputs already exist and this stage never overwrites:\n  "
                 + "\n  ".join(clashes)
                 + "\nBump window_scan.version_tag in workflow/config.yaml instead.")


def load_expected_snp_counts(path):
    """chrom -> n_snps plus the TOTAL row, from Stage 03's summary."""
    df = pd.read_csv(path, sep="\t")
    total = int(df.loc[df["chrom"] == "TOTAL", "n_snps"].iloc[0])
    per_chrom = {r.chrom: int(r.n_snps) for r in df.itertuples() if r.chrom != "TOTAL"}
    return per_chrom, total


def check_cohort(vcf_samples, samples_file):
    """Anti-silent-loss check at the Stage 03 -> Stage 04 boundary."""
    expected = [ln.strip() for ln in Path(samples_file).read_text().splitlines() if ln.strip()]
    if len(vcf_samples) != len(expected):
        sys.exit(f"ERROR: VCF carries {len(vcf_samples)} samples but {samples_file} "
                 f"lists {len(expected)}. The Stage 03 -> Stage 04 boundary has lost "
                 f"or gained samples; do not report this scan as the cohort scan.")
    if set(vcf_samples) != set(expected):
        only_vcf = sorted(set(vcf_samples) - set(expected))[:5]
        only_lst = sorted(set(expected) - set(vcf_samples))[:5]
        sys.exit(f"ERROR: same sample count but different sample sets.\n"
                 f"  only in VCF: {only_vcf}\n  only in list: {only_lst}")
    return len(expected)


# ---------------------------------------------------------------------------
# L4 — window de-duplication
# ---------------------------------------------------------------------------
def dedupe_windows(windows, pos_slices):
    """Collapse windows carrying an identical tuple of variant positions,
    keeping the FIRST occurrence (NB1 `np.unique(positions, return_index=True)`).

    windows     (n, 2) int array of [start, end]
    pos_slices  list of (lo, hi) index pairs into the chromosome's sorted `pos`;
                an identical index span is exactly an identical position tuple
    Returns the indices of the kept rows, in input (genomic) order.
    """
    seen = set()
    keep = []
    for i, span in enumerate(pos_slices):
        if span not in seen:
            seen.add(span)
            keep.append(i)
    return np.asarray(keep, dtype=int)


def _self_test_dedupe():
    """Windows 0 and 1 below hold the same two variants; window 2 holds one.
    Keeping the first occurrence must yield indices [0, 2]."""
    windows = np.array([[1, 200], [51, 250], [201, 400]])
    slices = [(0, 2), (0, 2), (2, 3)]
    got = dedupe_windows(windows, slices)
    assert got.tolist() == [0, 2], f"dedupe self-test failed: {got.tolist()}"


# ---------------------------------------------------------------------------
# Core-genome overlap (L3)
# ---------------------------------------------------------------------------
def load_core_coverage(bed_path, chrom_lengths):
    """chrom -> 1-based boolean coverage array over pk_core_coding.bed.

    BED is 0-based half-open; the returned array is indexed by 1-based position
    so cov[start:end+1] is exactly the window span.
    """
    cov = {c: np.zeros(L + 2, dtype=bool) for c, L in chrom_lengths.items()}
    with Path(bed_path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            chrom, s, e = f[0], int(f[1]), int(f[2])
            if chrom not in cov:
                continue
            lo, hi = s + 1, min(e, chrom_lengths[chrom])  # -> 1-based inclusive
            if hi >= lo:
                cov[chrom][lo:hi + 1] = True
    return cov


# ---------------------------------------------------------------------------
# Per-window statistics
# ---------------------------------------------------------------------------
def diversity(counts):
    """(het, entropy) from unique-haplotype sample counts. He = 1 - sum f^2;
    entropy = -sum f ln f — Siegel's definitions, unchanged."""
    total = counts.sum()
    if total == 0:
        return np.nan, np.nan
    f = counts / total
    return 1.0 - float(np.sum(f ** 2)), -float(np.sum(f * np.log(f)))


def window_stats(gt_block, hap_block):
    """All per-window statistics for one window.

    gt_block   (n_variants, n_samples, 2) int8, missing = -1. UN-haploidified (L2).
    hap_block  (n_variants, n_samples)    int8, missing = -1. Haploidified (L2).
    """
    # --- raw / Siegel-faithful branch (L1, L2) -----------------------------
    _, index, counts = np.unique(gt_block, axis=1, return_index=True, return_counts=True)
    het_raw, entropy = diversity(counts)
    n_alleles = int(len(counts))

    reps = gt_block[:, index, :]                       # (n_var, n_alleles, 2)
    with_missing = int(np.any(reps < 0, axis=(0, 2)).sum())
    a0, a1 = reps[:, :, 0], reps[:, :, 1]
    with_het = int(np.any((a0 != a1) & (a0 >= 0) & (a1 >= 0), axis=0).sum())

    # --- no-missing branch (L1) -------------------------------------------
    keep = ~np.any(hap_block < 0, axis=0)
    n_used = int(keep.sum())
    if n_used:
        _, counts_nm = np.unique(hap_block[:, keep], axis=1, return_counts=True)
        het_nm, entropy_nm = diversity(counts_nm)
        n_alleles_nm = int(len(counts_nm))
        delta = het_raw - het_nm
    else:
        het_nm = entropy_nm = delta = np.nan
        n_alleles_nm = 0

    return (n_alleles, het_raw, entropy, with_missing, with_het,
            n_alleles_nm, het_nm, entropy_nm, delta, n_used)


def scan_chromosome(chrom, pos, gt, cov, size, step):
    """Window rows for one chromosome, plus that chromosome's funnel counts."""
    n_variants, windows = allel.windowed_count(pos, size=size, step=step, start=1)

    non_empty = np.flatnonzero(n_variants != 0)          # NB1 keeps != 0, not > 1
    windows_ne = windows[non_empty]

    # Index span of the variants inside each surviving window. pos is sorted.
    slices = [(int(np.searchsorted(pos, s, "left")), int(np.searchsorted(pos, e, "right")))
              for s, e in windows_ne]

    keep = dedupe_windows(windows_ne, slices)            # L4
    windows_u = windows_ne[keep]
    slices_u = [slices[i] for i in keep]

    # One haploid realisation of the cohort per chromosome (L2). np.random is
    # seeded once by the caller and chromosomes are walked in a fixed order, so
    # the realisation is reproducible and is identical for every window that
    # overlaps a given variant.
    hap = allel.GenotypeArray(gt).haploidify_samples().values.astype(np.int8)

    rows = []
    for (ws, we), (lo, hi) in zip(windows_u, slices_u):
        stats = window_stats(gt[lo:hi], hap[lo:hi])
        span = cov[ws:we + 1]
        rows.append((
            chrom, int(ws), int(we), int((ws + we) // 2),
            hi - lo, ",".join(str(p) for p in pos[lo:hi]),
            *stats,
            float(span.mean()) if span.size else np.nan,
        ))

    funnel = {
        "windows_total": int(len(windows)),
        "windows_non_empty": int(len(windows_ne)),
        "windows_unique": int(len(windows_u)),
        "variants_read": int(len(pos)),
    }
    return rows, funnel


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    _self_test_dedupe()
    if args.missing_policy != "both":
        sys.exit(f"ERROR: --missing-policy={args.missing_policy!r} is not implemented. "
                 "L1 fixes this to 'both' (emit het_raw and het_nomiss, filter on "
                 "het_raw). Change the lock in handoffs/04_window_scan.md before "
                 "changing the config key.")
    refuse_to_overwrite([args.out_windows, args.out_heterozygome, args.out_summary])

    print(f"==> window_scan: size={args.window_size} step={args.step} "
          f"min_snps={args.min_snps} min_het={args.min_het} seed={args.seed}")

    print(f"    reading {args.vcf}")
    callset = allel.read_vcf(
        args.vcf, fields=["variants/CHROM", "variants/POS", "calldata/GT", "samples"])
    if callset is None:
        sys.exit(f"ERROR: no records read from {args.vcf}")

    chroms = callset["variants/CHROM"]
    pos_all = callset["variants/POS"]
    gt_all = callset["calldata/GT"]
    n_samples = check_cohort(list(callset["samples"]), args.samples)
    print(f"    cohort: {n_samples} samples (matches {args.samples})")

    expected_per_chrom, expected_total = load_expected_snp_counts(args.snps_summary)
    if len(pos_all) != expected_total:
        sys.exit(f"ERROR: read {len(pos_all)} variants but {args.snps_summary} "
                 f"says {expected_total}.")
    print(f"    variants: {len(pos_all):,} (matches {args.snps_summary} TOTAL)")

    observed = {c: int((chroms == c).sum()) for c in np.unique(chroms)}
    if observed != expected_per_chrom:
        diff = {c: (observed.get(c), expected_per_chrom.get(c))
                for c in set(observed) | set(expected_per_chrom)
                if observed.get(c) != expected_per_chrom.get(c)}
        sys.exit(f"ERROR: per-chromosome variant counts disagree with "
                 f"{args.snps_summary} (observed, expected): {diff}")
    order = sorted(expected_per_chrom)          # deterministic; L2 depends on it
    print(f"    per-chromosome counts match across {len(order)} chromosomes")

    # frac_in_core needs a 1-based coverage array per chromosome (L3). Length is
    # the furthest of the last variant's window and the last BED interval.
    bed_ends = {}
    with Path(args.bed).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            bed_ends[f[0]] = max(bed_ends.get(f[0], 0), int(f[2]))
    chrom_lengths = {c: max(int(pos_all[chroms == c].max()) + args.window_size,
                            bed_ends.get(c, 0))
                     for c in order}
    cov = load_core_coverage(args.bed, chrom_lengths)
    print(f"    core coverage loaded from {args.bed}")

    np.random.seed(args.seed)                   # L2 — seeded once, then `order`
    rows, per_chrom_funnel = [], {}
    for c in order:
        m = chroms == c
        r, f = scan_chromosome(c, pos_all[m], gt_all[m], cov[c],
                               args.window_size, args.step)
        rows.extend(r)
        per_chrom_funnel[c] = f
        print(f"    {c}: {f['variants_read']:,} SNPs -> {f['windows_non_empty']:,} "
              f"non-empty -> {f['windows_unique']:,} unique windows")

    df = pd.DataFrame(rows, columns=COLUMNS)

    # --- funnel ------------------------------------------------------------
    tot = {k: sum(f[k] for f in per_chrom_funnel.values())
           for k in ("windows_total", "windows_non_empty", "windows_unique")}
    n_min_snps = int((df["n_variants"] >= args.min_snps).sum())
    het_mask = (df["n_variants"] >= args.min_snps) & (df["het_raw"] >= args.min_het)
    n_candidates = int(het_mask.sum())
    preview = int((df["n_variants"].between(3, 10) & (df["het_raw"] >= 0.60)).sum())

    Path(args.out_windows).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_windows, sep="\t", index=False, compression="gzip",
              float_format="%.6g")
    df[het_mask].to_csv(args.out_heterozygome, sep="\t", index=False,
                        float_format="%.6g")

    summary = [("scope", "metric", "value")]
    for k, v in (("vcf", args.vcf), ("samples_file", args.samples),
                 ("bed", args.bed), ("n_samples", n_samples),
                 ("window_size_bp", args.window_size), ("step_bp", args.step),
                 ("min_snps", args.min_snps), ("min_heterozygosity", args.min_het),
                 ("haploidify_seed", args.seed), ("missing_policy", args.missing_policy),
                 ("scikit_allel_version", allel.__version__)):
        summary.append(("params", k, str(v)))
    for k, v in (("variants_read", len(pos_all)),
                 ("windows_total", tot["windows_total"]),
                 ("windows_non_empty", tot["windows_non_empty"]),
                 ("windows_unique_post_dedup", tot["windows_unique"]),
                 ("windows_ge_min_snps", n_min_snps),
                 ("heterozygome_ge_min_snps_and_min_het", n_candidates),
                 ("preview_3to10_snps_het_ge_0.60", preview)):
        summary.append(("funnel", k, str(v)))
    for c in order:
        f = per_chrom_funnel[c]
        for k in ("variants_read", "windows_non_empty", "windows_unique"):
            summary.append((c, k, str(f[k])))
        sub = df[df["chrom"] == c]
        summary.append((c, "heterozygome_windows",
                        str(int(((sub["n_variants"] >= args.min_snps)
                                 & (sub["het_raw"] >= args.min_het)).sum()))))
    for k, v in (("het_raw_median", df["het_raw"].median()),
                 ("het_nomiss_median", df["het_nomiss"].median()),
                 ("delta_het_median", df["delta_het"].median()),
                 ("frac_in_core_median", df["frac_in_core"].median()),
                 ("frac_in_core_eq_1_fraction", float((df["frac_in_core"] >= 1.0).mean()))):
        summary.append(("distribution", k, f"{v:.6g}"))

    with Path(args.out_summary).open("w") as fh:
        for r in summary:
            fh.write("\t".join(r) + "\n")

    print(f"\n==> funnel: {tot['windows_total']:,} total -> "
          f"{tot['windows_non_empty']:,} non-empty -> "
          f"{tot['windows_unique']:,} unique -> "
          f"{n_min_snps:,} >={args.min_snps} SNPs -> "
          f"{n_candidates:,} heterozygome ({preview:,} at 3-10 SNPs & He>=0.60)")
    print(f"==> delta_het median {df['delta_het'].median():.4f} "
          f"(het_raw {df['het_raw'].median():.4f} vs "
          f"het_nomiss {df['het_nomiss'].median():.4f})")
    print(f"==> wrote {args.out_windows}\n            {args.out_heterozygome}"
          f"\n            {args.out_summary}")


if __name__ == "__main__":
    main()
