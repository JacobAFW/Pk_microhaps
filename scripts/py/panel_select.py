#!/usr/bin/env python
"""
panel_select.py — Stage 05: choose a microhaplotype panel from the heterozygome.

Port of Siegel et al's notebook 2 selection step. Takes the Stage 04
heterozygome, applies the tighter selection filters, collapses overlapping
windows into independent loci, allocates markers across chromosomes by length,
and writes a panel of the requested size by one of two algorithms.

**STATUS: authored and linted, NEVER RUN.** Stage 04's outputs exist, so this
will run, but no result from it has been produced or reviewed. Treat the exit
criteria in `workflow/rules/05_panel_select.README.md` as untested.

---------------------------------------------------------------------------
THE ONE THING THAT IS EASY TO GET WRONG HERE
---------------------------------------------------------------------------
Stage 04 windows are 200 bp wide at a 50 bp step, so neighbouring candidate
windows overlap by up to 75% and describe THE SAME PIECE OF GENOME. Selecting
over raw windows would put four near-identical markers on one locus and call it
four markers.

Measured on the Stage 04 output (`runs/2026-09-25_stage04-window-scan/`,
phase-4 review):

    selection-eligible windows      19,674
    non-overlapping merged loci      5,915      (3.33 windows per locus)

So the real headroom for a 100-marker panel is ~59x, not ~197x. This script
therefore collapses overlapping windows to loci FIRST, keeps the best window per
locus, and then enforces `panel_select.min_spacing_bp` between chosen markers so
two adjacent loci cannot both be taken. Both counts are reported.

---------------------------------------------------------------------------
Usage:
    python scripts/py/panel_select.py \
        --heterozygome outputs/scan/v1_2026-09-25/heterozygome.tsv \
        --fai data/reference/strain_preDB_version/strain_A1_H.1.Icor.fasta.fai \
        --min-snps 3 --max-snps 10 --min-het 0.60 \
        --min-spacing-bp 10000 \
        --algorithm greedy --panel-size 100 \
        --out-candidates outputs/panel/v1/candidates.tsv \
        --out-panel outputs/panel/v1/panel_greedy_100.tsv

Algorithms:
    greedy          highest het_raw first, subject to spacing and the
                    per-chromosome allocation.
    evenly_spaced   allocate per chromosome by length, place ideal positions at
                    equal intervals, take the nearest eligible locus to each.

Exit codes:
    0  panel written
    1  a guard failed (message on stderr says which)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PANEL_COLUMNS = [
    "chrom", "start", "end", "midpoint", "n_snps", "het_raw", "het_nomiss",
    "entropy", "effective_cardinality", "frac_in_core", "variant_positions",
    "locus_id", "rank_in_chrom",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--heterozygome", required=True, help="Stage 04 heterozygome.tsv")
    p.add_argument("--fai", required=True,
                   help=".fai whose contig names match the VCF (ordered_PKNH_*_v2), "
                        "used for length-proportional per-chromosome allocation")
    p.add_argument("--min-snps", type=int, required=True)
    p.add_argument("--max-snps", type=int, required=True)
    p.add_argument("--min-het", type=float, required=True)
    p.add_argument("--min-spacing-bp", type=int, required=True,
                   help="minimum distance between selected markers; also stops two "
                        "overlapping windows being taken as separate markers")
    p.add_argument("--algorithm", required=True, choices=["greedy", "evenly_spaced"])
    p.add_argument("--panel-size", type=int, required=True)
    p.add_argument("--out-candidates", required=True)
    p.add_argument("--out-panel", required=True)
    return p.parse_args()


def refuse_to_overwrite(paths):
    """Standing rule: nothing existing gets overwritten."""
    clashes = [p for p in paths if Path(p).exists()]
    if clashes:
        sys.exit("ERROR: these outputs already exist and this stage never "
                 "overwrites:\n  " + "\n  ".join(clashes)
                 + "\nBump panel_select.version_tag in workflow/config.yaml instead.")


def load_contig_lengths(fai_path):
    """chrom -> length, nuclear contigs only (the scan covers no others)."""
    lengths = {}
    with Path(fai_path).open() as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                lengths[f[0]] = int(f[1])
    return lengths


def merge_to_loci(df):
    """Collapse overlapping/bookending windows into independent loci.

    Stage 04's 200 bp / 50 bp tiling makes neighbouring candidate windows
    describe the same genome. Returns `df` with a `locus_id` column; loci are
    numbered per chromosome in genomic order.
    """
    out = []
    for chrom, sub in df.groupby("chrom", sort=True):
        sub = sub.sort_values("window_start").copy()
        ids, cur, end = [], -1, -1
        for s, e in zip(sub.window_start, sub.window_end):
            if s > end:
                cur += 1
                end = e
            else:
                end = max(end, e)
            ids.append(f"{chrom}:L{cur:05d}")
        sub["locus_id"] = ids
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def best_per_locus(df):
    """One representative window per locus: highest het_raw, then most SNPs,
    then lowest start — a total order, so the choice is deterministic."""
    return (df.sort_values(["locus_id", "het_raw", "n_variants", "window_start"],
                           ascending=[True, False, False, True])
              .groupby("locus_id", as_index=False)
              .first())


def allocate_by_length(lengths, chroms, panel_size):
    """Largest-remainder allocation of `panel_size` markers across chromosomes,
    proportional to length. Guarantees the parts sum to the whole."""
    lens = {c: lengths[c] for c in chroms if c in lengths}
    if not lens:
        sys.exit("ERROR: no chromosome in the heterozygome matched the .fai. "
                 "Wrong .fai? The scan uses ordered_PKNH_*_v2 naming.")
    total = sum(lens.values())
    exact = {c: panel_size * L / total for c, L in lens.items()}
    alloc = {c: int(np.floor(v)) for c, v in exact.items()}
    remainder = panel_size - sum(alloc.values())
    for c, _ in sorted(exact.items(), key=lambda kv: kv[1] - np.floor(kv[1]),
                       reverse=True)[:remainder]:
        alloc[c] += 1
    return alloc


def pick_greedy(cand, n, min_spacing):
    """Highest het_raw first; skip anything within `min_spacing` of a pick."""
    chosen = []
    for r in cand.sort_values(["het_raw", "n_variants"],
                              ascending=[False, False]).itertuples():
        if len(chosen) >= n:
            break
        if all(abs(r.midpoint - c.midpoint) >= min_spacing for c in chosen):
            chosen.append(r)
    return chosen


def pick_evenly_spaced(cand, n, min_spacing, chrom_length):
    """Place n ideal positions at equal intervals and take the nearest eligible
    locus to each, skipping any that would violate the spacing rule."""
    if n <= 0 or cand.empty:
        return []
    targets = np.linspace(chrom_length / (2 * n), chrom_length - chrom_length / (2 * n), n)
    chosen, used = [], set()
    for t in targets:
        order = (cand.midpoint - t).abs().sort_values().index
        for idx in order:
            if idx in used:
                continue
            r = cand.loc[idx]
            if all(abs(r.midpoint - c.midpoint) >= min_spacing for c in chosen):
                chosen.append(r)
                used.add(idx)
                break
    return chosen


def main():
    args = parse_args()
    refuse_to_overwrite([args.out_candidates, args.out_panel])

    het = pd.read_csv(args.heterozygome, sep="\t")
    if het.empty:
        sys.exit(f"ERROR: {args.heterozygome} is empty — Stage 04 produced no "
                 "candidate windows, so there is nothing to select from.")
    print(f"==> {len(het):,} heterozygome windows from {args.heterozygome}")

    # --- selection filters (tighter than the scan's) ----------------------
    sel = het[(het.n_variants >= args.min_snps)
              & (het.n_variants <= args.max_snps)
              & (het.het_raw >= args.min_het)].copy()
    print(f"    {len(sel):,} pass selection filters "
          f"({args.min_snps}-{args.max_snps} SNPs, het_raw >= {args.min_het})")
    if sel.empty:
        sys.exit("ERROR: no window passes the selection filters.")

    # --- collapse overlapping windows to independent loci ------------------
    sel = merge_to_loci(sel)
    cand = best_per_locus(sel)
    n_loci = cand.locus_id.nunique()
    print(f"    {len(sel):,} windows -> {n_loci:,} non-overlapping loci "
          f"({len(sel) / n_loci:.2f} windows per locus)")
    print(f"    headroom for a {args.panel_size}-marker panel: "
          f"{n_loci / args.panel_size:.0f}x")
    if n_loci < args.panel_size:
        sys.exit(f"ERROR: only {n_loci} independent loci available but "
                 f"{args.panel_size} markers requested.")

    cand = cand.rename(columns={"window_start": "start", "window_end": "end",
                                "n_variants": "n_snps"})
    # Effective number of haplotypes: 1 / sum(f^2) = 1 / (1 - He).
    cand["effective_cardinality"] = 1.0 / (1.0 - cand.het_raw).clip(lower=1e-12)

    Path(args.out_candidates).parent.mkdir(parents=True, exist_ok=True)
    cand.sort_values(["chrom", "start"]).to_csv(
        args.out_candidates, sep="\t", index=False, float_format="%.6g")
    print(f"==> wrote {args.out_candidates}")

    # --- per-chromosome allocation, proportional to length ----------------
    lengths = load_contig_lengths(args.fai)
    chroms = sorted(cand.chrom.unique())
    alloc = allocate_by_length(lengths, chroms, args.panel_size)
    print(f"    allocation across {len(alloc)} chromosomes: "
          + ", ".join(f"{c.replace('ordered_PKNH_', '').replace('_v2', '')}={n}"
                      for c, n in sorted(alloc.items())))

    # --- select --------------------------------------------------------------
    picks, shortfall = [], 0
    for chrom in chroms:
        want = alloc.get(chrom, 0)
        if want == 0:
            continue
        sub = cand[cand.chrom == chrom]
        if args.algorithm == "greedy":
            got = pick_greedy(sub, want, args.min_spacing_bp)
        else:
            got = pick_evenly_spaced(sub, want, args.min_spacing_bp, lengths[chrom])
        if len(got) < want:
            shortfall += want - len(got)
            print(f"    WARNING: {chrom} wanted {want}, got {len(got)} "
                  f"(spacing {args.min_spacing_bp:,} bp is binding)")
        for rank, r in enumerate(got, 1):
            row = r._asdict() if hasattr(r, "_asdict") else dict(r)
            row["rank_in_chrom"] = rank
            picks.append(row)

    panel = pd.DataFrame(picks)
    if panel.empty:
        sys.exit("ERROR: no markers selected.")
    for c in PANEL_COLUMNS:
        if c not in panel.columns:
            panel[c] = np.nan
    panel = panel[PANEL_COLUMNS].sort_values(["chrom", "start"])
    panel.to_csv(args.out_panel, sep="\t", index=False, float_format="%.6g")

    print(f"\n==> panel: {len(panel)} markers "
          f"({args.algorithm}, requested {args.panel_size})")
    if shortfall:
        print(f"    SHORTFALL {shortfall} — relax panel_select.min_spacing_bp "
              f"or reduce the panel size.")
    print(f"    het_raw  median {panel.het_raw.median():.4f} "
          f"[{panel.het_raw.min():.4f}, {panel.het_raw.max():.4f}]")
    print(f"    n_snps   median {panel.n_snps.median():.0f}")
    gaps = (panel.sort_values(['chrom', 'start'])
                 .groupby("chrom").midpoint.diff().dropna())
    if len(gaps):
        print(f"    spacing  median {gaps.median():,.0f} bp, min {gaps.min():,.0f} bp")
    print(f"==> wrote {args.out_panel}")


if __name__ == "__main__":
    main()
