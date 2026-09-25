#!/usr/bin/env python3
"""make_smiss.py — write a PLINK2-shaped .smiss/.fam pair for find_duplicates.R.

Stage 02b, rule remove_duplicates. The vendored find_duplicates.R is kept
byte-identical to the agnostic upstream, so it expects exactly what PLINK2
produces: a `.fam` (headerless, FID IID PAT MAT SEX PHENO) and a `.smiss`
(header `#FID  IID  MISSING_CT  OBS_CT  F_MISS`).

This repo has no plink2 — Stages 01-03 are pure bcftools — so we synthesise
both from Stage 02's already-computed per-sample callable fraction:

    F_MISS = 1 - callable_fraction

That is missingness over pk_core_coding.bed, whereas agnostic measured it over
its `cleaned` bfile. The number differs; the *ordering* is what find_duplicates.R
uses, and it only ever compares two technical replicates of the same sample. One
tie-break convention, no extra pass over a 17 GiB VCF.

FID is set equal to IID, matching agnostic's `plink2 --double-id`, so the
`--remove` list find_duplicates.R writes is in "sample sample" form.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--callable", required=True, type=Path,
                    help="outputs/qc/callable_fraction.tsv from Stage 02")
    ap.add_argument("--samples", required=True, type=Path,
                    help="cohort sample list to restrict to, one ID per line")
    ap.add_argument("--out-fam", required=True, type=Path)
    ap.add_argument("--out-smiss", required=True, type=Path)
    args = ap.parse_args()

    wanted = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    if not wanted:
        sys.exit(f"ERROR: {args.samples} is empty")

    stats: dict[str, tuple[int, int, float]] = {}
    with args.callable.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_called = int(row["n_called"])
            n_total = int(row["n_total"])
            f_miss = 1.0 - float(row["callable_fraction"])
            stats[row["sample"]] = (n_total - n_called, n_total, f_miss)

    missing = [s for s in wanted if s not in stats]
    if missing:
        sys.exit(f"ERROR: {len(missing)} cohort sample(s) absent from "
                 f"{args.callable}, e.g. {missing[:5]}")

    args.out_fam.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fam.open("w") as fh:
        for s in wanted:
            fh.write(f"{s}\t{s}\t0\t0\t0\t-9\n")

    with args.out_smiss.open("w") as fh:
        fh.write("#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\n")
        for s in wanted:
            miss_ct, obs_ct, f_miss = stats[s]
            fh.write(f"{s}\t{s}\t{miss_ct}\t{obs_ct}\t{f_miss:.6f}\n")

    f_miss_vals = [stats[s][2] for s in wanted]
    print(f"[make_smiss] samples: {len(wanted)}")
    print(f"[make_smiss] F_MISS: min {min(f_miss_vals):.4f}  "
          f"median {sorted(f_miss_vals)[len(f_miss_vals) // 2]:.4f}  "
          f"max {max(f_miss_vals):.4f}")
    print(f"[make_smiss] wrote {args.out_fam} and {args.out_smiss}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
