#!/usr/bin/env python3
"""compare_fws_sources.py — the A1 verdict.

Deviation A1 (CONTEXT-method-deviations.md, open since 2026-05-12): our Fws is
imported from the Pop-gen pipeline's `Proportion` column and has never been
verified. moimix was not installed, so there was no cross-check and the decision
sat open.

This compares the imported `Proportion` against Fws recomputed by
scripts/R/run_moimix.R on the same cohort, and reports three numbers:

  1. Pearson r
  2. max absolute difference
  3. samples that cross the Fws threshold under one source but not the other

Number 3 IS the verdict. If it is zero, `Proportion` is fit for purpose whatever
it is called and A1 closes — the discovery cohort is unaffected by the choice of
source. If it is non-zero, A1 is a real problem: the cohort membership itself
depends on which column you believe, and the exit code says so.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def read_fws(path: Path, sample_col: str | None, value_col: str | None) -> dict[str, float]:
    """Read a two-column-ish TSV of per-sample Fws, guessing columns if needed."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields:
            sys.exit(f"ERROR: no header in {path}")
        scol = sample_col or next(
            (c for c in ("sample", "Sample", "sample_id", "IID") if c in fields), fields[0])
        vcol = value_col or next(
            (c for c in ("Proportion", "fws", "Fws", "FWS") if c in fields), fields[-1])
        out: dict[str, float] = {}
        for row in reader:
            raw = (row.get(vcol) or "").strip()
            if raw in ("", "NA", "NaN"):
                continue
            out[(row.get(scol) or "").strip()] = float(raw)
    print(f"[compare_fws] {path.name}: {len(out)} samples "
          f"(sample col {scol!r}, value col {vcol!r})")
    return out


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imported", required=True, type=Path,
                    help="Pop-gen fws table (the `Proportion` column)")
    ap.add_argument("--computed", required=True, type=Path,
                    help="fws recomputed here by run_moimix.R")
    ap.add_argument("--threshold", type=float, default=0.95)
    ap.add_argument("--out-report", required=True, type=Path)
    ap.add_argument("--fail-on-crossing", action="store_true",
                    help="exit non-zero if any sample crosses the threshold "
                         "differently between sources (handoff #05 stop condition)")
    args = ap.parse_args()

    imported = read_fws(args.imported, None, None)
    computed = read_fws(args.computed, None, None)

    shared = sorted(set(imported) & set(computed))
    if not shared:
        sys.exit("ERROR: no overlapping sample IDs between the two Fws sources")

    xs = [imported[s] for s in shared]
    ys = [computed[s] for s in shared]
    r = pearson(xs, ys)
    diffs = [abs(a - b) for a, b in zip(xs, ys)]
    max_abs = max(diffs)

    crossings = [
        (s, imported[s], computed[s])
        for s in shared
        if (imported[s] >= args.threshold) != (computed[s] >= args.threshold)
    ]

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"n_shared_samples\t{len(shared)}\n")
        fh.write(f"n_imported_only\t{len(set(imported) - set(computed))}\n")
        fh.write(f"n_computed_only\t{len(set(computed) - set(imported))}\n")
        fh.write(f"pearson_r\t{r:.6f}\n")
        fh.write(f"max_abs_difference\t{max_abs:.6f}\n")
        fh.write(f"mean_abs_difference\t{sum(diffs) / len(diffs):.6f}\n")
        fh.write(f"threshold\t{args.threshold}\n")
        fh.write(f"n_threshold_crossings\t{len(crossings)}\n")
        fh.write(f"a1_verdict\t{'CLOSED' if not crossings else 'OPEN'}\n")
        if crossings:
            fh.write("\nsample\timported\tcomputed\n")
            for s, a, b in crossings:
                fh.write(f"{s}\t{a:.6f}\t{b:.6f}\n")

    print(f"\n=== A1 verdict ===")
    print(f"shared samples:        {len(shared)}")
    print(f"Pearson r:             {r:.6f}")
    print(f"max abs difference:    {max_abs:.6f}")
    print(f"threshold crossings:   {len(crossings)}   <-- THE VERDICT")
    for s, a, b in crossings:
        print(f"    {s}: imported={a:.6f} computed={b:.6f}")
    print(f"verdict:               {'CLOSED' if not crossings else 'OPEN'}")
    print(f"wrote {args.out_report}")

    if crossings and args.fail_on_crossing:
        sys.exit(f"ERROR: {len(crossings)} sample(s) cross the Fws {args.threshold} "
                 f"threshold differently between sources. A1 is a real problem, not "
                 f"a formality — it invalidates the discovery cohort. Stop and flag.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
