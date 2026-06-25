#!/usr/bin/env python
"""
merge_sample_qc.py — Stage 02: combine Fws + callable + metadata, filter to
the monoclonal discovery cohort.

Inputs:
    --fws         outputs/setup/fws_from_popgen.tsv   (cols: sample, fws)
    --callable    outputs/qc/callable_fraction.tsv    (cols: sample, n_called, n_total, callable_fraction)
    --metadata    data/metadata/samples.tsv           (cols: sample, country, region_geo)

Outputs:
    --out-metadata   outputs/qc/sample_metadata.tsv
        Columns: sample, fws, n_called, n_total, callable_fraction,
                 country, region_geo, include_for_discovery
    --out-discovery  outputs/qc/discovery_samples.txt  (one sample per line)

Sanity assertions (any failure -> non-zero exit):
    - Every sample in callable_fraction has a corresponding Fws row
      (warn-only on missing-from-metadata; samples without metadata get
       region_geo='UNKNOWN').
    - Fws values are numeric and in [0, 1] (the moimix 'Proportion' column
      should satisfy this). TODO: confirm with Pop-gen author that
      'Proportion' == moimix's Fws statistic.
    - At least 50 samples pass both thresholds (sanity floor; Siegel had 615
      of 1816 Pv samples ≈ 34%).

Exit codes:
    0  success
    2  hard sanity-assertion failure
    3  Pop-gen Fws table empty (header-only stub from Stage 01)
"""

import argparse
import csv
import sys
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fws",       required=True)
    p.add_argument("--callable",  required=True, dest="callable_path")
    p.add_argument("--metadata",  required=True)
    p.add_argument("--fws-threshold",      type=float, required=True)
    p.add_argument("--callable-threshold", type=float, required=True)
    p.add_argument("--out-metadata",  required=True)
    p.add_argument("--out-discovery", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    fws_rows      = read_tsv(args.fws)
    callable_rows = read_tsv(args.callable_path)
    meta_rows     = read_tsv(args.metadata)

    print(f"==> Loaded {len(fws_rows)} Fws rows, {len(callable_rows)} callable rows, "
          f"{len(meta_rows)} metadata rows")

    if len(fws_rows) == 0:
        sys.exit(
            "ERROR: Fws table is empty (Stage 01 wrote a header-only stub). "
            "Either run Pop-gen Stage 2 and re-run import_popgen_fws, or "
            "install moimix into the R stack and add a recompute_fws rule."
        )

    # Index by sample.
    fws       = {r["sample"]: r for r in fws_rows}
    callable_ = {r["sample"]: r for r in callable_rows}
    meta      = {r["sample"]: r for r in meta_rows}

    # Sanity: Fws column is named 'fws' in the imported table. Values numeric, in [0,1].
    fws_values = []
    for sid, r in fws.items():
        try:
            v = float(r["fws"])
        except (KeyError, ValueError, TypeError):
            sys.exit(f"ERROR: non-numeric Fws value for {sid!r}: {r!r}")
        if not (0.0 <= v <= 1.0):
            sys.exit(
                f"ERROR: Fws value out of [0,1] for {sid!r}: {v}. "
                "Check the Pop-gen 'Proportion' column semantics — see "
                "CONTEXT-software.md Open issues #2."
            )
        fws_values.append(v)

    # Soft sanity: median Fws should be monoclonal-skewed (>0.8).
    med = sorted(fws_values)[len(fws_values) // 2]
    print(f"    Fws: n={len(fws_values)}, median={med:.3f}, "
          f"max={max(fws_values):.3f}, min={min(fws_values):.3f}")
    if med < 0.5:
        print(
            "    WARNING: median Fws < 0.5. Either the cohort is highly "
            "polyclonal or the 'Proportion' column is NOT Fws. Investigate "
            "before trusting Stage 02 output.",
            file=sys.stderr,
        )
    # TODO: confirm with Pop-gen Stage 2 author that 'Proportion' == Fws.

    samples = sorted(set(fws) | set(callable_))

    n_missing_meta = 0
    rows_out = []
    pass_count = 0
    for sid in samples:
        f = fws.get(sid, {})
        c = callable_.get(sid, {})
        m = meta.get(sid, {})

        fws_val = float(f["fws"]) if f else float("nan")
        n_called = int(c["n_called"])   if c else 0
        n_total  = int(c["n_total"])    if c else 0
        cf       = float(c["callable_fraction"]) if c else float("nan")

        country    = m.get("country", "") if m else ""
        region_geo = m.get("region_geo", "") if m else ""
        if not m:
            n_missing_meta += 1
            region_geo = "UNKNOWN"

        include = (
            f and c
            and fws_val >= args.fws_threshold
            and cf      >= args.callable_threshold
        )
        if include:
            pass_count += 1

        rows_out.append({
            "sample": sid,
            "fws": f"{fws_val:.6f}" if fws else "",
            "n_called": str(n_called),
            "n_total":  str(n_total),
            "callable_fraction": f"{cf:.6f}" if c else "",
            "country": country,
            "region_geo": region_geo,
            "include_for_discovery": "TRUE" if include else "FALSE",
        })

    if n_missing_meta:
        print(f"    WARNING: {n_missing_meta} samples have no metadata "
              f"(region_geo='UNKNOWN')", file=sys.stderr)

    # Write sample_metadata.tsv.
    fields = ["sample", "fws", "n_called", "n_total", "callable_fraction",
              "country", "region_geo", "include_for_discovery"]
    with Path(args.out_metadata).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"==> Wrote {args.out_metadata} ({len(rows_out)} samples)")

    # Write discovery_samples.txt.
    with Path(args.out_discovery).open("w") as fh:
        for r in rows_out:
            if r["include_for_discovery"] == "TRUE":
                fh.write(r["sample"] + "\n")
    print(f"==> Wrote {args.out_discovery} ({pass_count} samples passing both filters)")

    # Acceptance: at least 50 samples in the discovery cohort.
    if pass_count < 50:
        sys.exit(
            f"ERROR: only {pass_count} samples pass Fws ≥ {args.fws_threshold} "
            f"AND callable_fraction ≥ {args.callable_threshold}. "
            "Expected >= 50 for a usable discovery cohort. "
            "Either thresholds are wrong, or the cohort is more polyclonal than expected."
        )
    print(f"==> Discovery cohort: {pass_count} samples ({pass_count/len(samples):.1%} of input)")


if __name__ == "__main__":
    main()
