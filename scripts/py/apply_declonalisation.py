#!/usr/bin/env python3
"""apply_declonalisation.py — thin each clonal group to one representative.

Stage 02b, rule declonalised_cohort (decision L3).

He = 1 - sum(f^2), so a repeated haplotype pushes sum(f^2) up and He down. The
direction is conservative — you miss good windows rather than selecting bad ones
— but the deflation concentrates in windows where the repeated genotype is
distinctive, so it distorts the *ranking*. Panel selection is a ranking exercise
(top-N per region), so clonal over-representation yields a suboptimal panel, not
a wrong one. It also bites upstream at Stage 03's MAF >= 0.10 gate, changing
which SNPs enter at all.

With --declonalize the cohort keeps one representative per clonal group, chosen
by lowest missingness — the same tie-break find_duplicates.R uses for technical
replicates. One convention, not two.

Without it this is a pass-through copy, so the DAG shape is identical either way
and the sensitivity analysis is one config key away.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def read_f_miss(path: Path) -> dict[str, float]:
    f_miss: dict[str, float] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            iid = row.get("IID")
            if iid is not None:
                f_miss[iid] = float(row["F_MISS"])
    return f_miss


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, type=Path,
                    help="post-dedupe cohort, one ID per line")
    ap.add_argument("--clonal-clusters", required=True, type=Path,
                    help="outputs/qc/clonal_clusters.tsv from clonal_clusters.R")
    ap.add_argument("--smiss", required=True, type=Path,
                    help="per-sample missingness, for the representative tie-break")
    ap.add_argument("--declonalize", action="store_true")
    ap.add_argument("--out-samples", required=True, type=Path)
    ap.add_argument("--out-report", required=True, type=Path)
    args = ap.parse_args()

    cohort = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    cohort_set = set(cohort)
    f_miss = read_f_miss(args.smiss)

    groups: dict[str, list[str]] = defaultdict(list)
    with args.clonal_clusters.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            grp = (row.get("clonal_group") or "").strip()
            sample = (row.get("sample") or "").strip()
            if grp and grp != "NA" and sample in cohort_set:
                groups[grp].append(sample)

    n_grouped = sum(len(v) for v in groups.values())
    print(f"[declonalise] cohort in:        {len(cohort)}")
    print(f"[declonalise] clonal groups:    {len(groups)}")
    print(f"[declonalise] samples in group: {n_grouped}")

    dropped: list[tuple[str, str, str]] = []
    if args.declonalize:
        for grp, members in sorted(groups.items()):
            # Lowest missingness wins; ID as a deterministic tie-break so the
            # representative does not depend on dict ordering.
            ranked = sorted(members, key=lambda s: (f_miss.get(s, 1.0), s))
            keep = ranked[0]
            for loser in ranked[1:]:
                dropped.append((loser, grp, keep))
        drop_set = {d[0] for d in dropped}
        kept = [s for s in cohort if s not in drop_set]
    else:
        kept = list(cohort)

    args.out_samples.parent.mkdir(parents=True, exist_ok=True)
    args.out_samples.write_text("".join(f"{s}\n" for s in kept))

    with args.out_report.open("w") as fh:
        fh.write("sample\tclonal_group\trepresentative\tf_miss\taction\n")
        for grp, members in sorted(groups.items()):
            ranked = sorted(members, key=lambda s: (f_miss.get(s, 1.0), s))
            rep = ranked[0]
            for s in ranked:
                if not args.declonalize:
                    action = "kept (declonalize=false)"
                elif s == rep:
                    action = "kept (representative)"
                else:
                    action = "dropped"
                fh.write(f"{s}\t{grp}\t{rep}\t{f_miss.get(s, float('nan')):.6f}\t"
                         f"{action}\n")

    pct = 100.0 * len(dropped) / len(cohort) if cohort else 0.0
    print(f"[declonalise] declonalize:      {args.declonalize}")
    print(f"[declonalise] dropped:          {len(dropped)} ({pct:.1f}% of cohort)")
    print(f"[declonalise] cohort out:       {len(kept)}")
    print(f"[declonalise] wrote {args.out_samples}")

    if args.declonalize and pct > 20.0:
        # Stop condition from handoff #05: more clonality than the discovery
        # design assumes, and a finding worth a conversation before it silently
        # reshapes the panel.
        sys.exit(f"ERROR: declonalisation dropped {pct:.1f}% of the cohort "
                 f"(> 20% stop condition). Sabah Pk may be more clonal than the "
                 f"discovery design assumes — inspect before proceeding.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
