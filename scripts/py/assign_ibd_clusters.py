#!/usr/bin/env python3
"""assign_ibd_clusters.py — map the discovery cohort onto imported ADMIXTURE clusters.

Stage 02b, decision L2. hmmIBD uses population allele frequencies, so IBD must be
computed within genetically homogeneous groups; pooling Indonesia + Sabah inflates
apparent IBD between samples from the same region (population structure
masquerading as relatedness).

We import cluster *labels* from the agnostic pop-gen pipeline
(data/metadata/admix_clusters.tsv), never a pairwise IBD result. A label is a
stable per-sample annotation; importing one keeps the regression comparison
apples-to-apples without importing the answer.

Joins are exact-first, then normalised (the same normalisation find_duplicates.R
applies: strip `duplicate_id_pattern`, then drop underscores and hyphens).
Samples with no label land in a residual group which is written to the report but
only gets a keep-list if --residual-clonal-calls is passed.

Outputs
  --out-assignments  sample <TAB> cluster   (every input sample, one row each)
  --out-summary      cluster <TAB> n_samples <TAB> ibd_eligible <TAB> reason
  --out-keep-dir     one <cluster>.keep per IBD-eligible cluster, "FID IID" format
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path


def normalise(sample: str, pattern: str) -> str:
    """Match find_duplicates.R: str_remove(pattern) then drop _ and -."""
    out = re.sub(pattern, "", sample) if pattern else sample
    return out.replace("_", "").replace("-", "")


def read_cluster_table(path: Path) -> "OrderedDict[str, str]":
    """Read admix_clusters.tsv -> {sample: cluster}, tolerating column aliases."""
    labels: "OrderedDict[str, str]" = OrderedDict()
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sample = row.get("Sample") or row.get("sample_id") or row.get("sample")
            cluster = row.get("Cluster") or row.get("cluster")
            if sample and cluster:
                labels[sample] = cluster
    if not labels:
        sys.exit(f"ERROR: no (Sample, Cluster) rows parsed from {path}")
    return labels


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, type=Path,
                    help="cohort sample list, one ID per line")
    ap.add_argument("--cluster-table", required=True, type=Path)
    ap.add_argument("--id-pattern", default="_DK.*",
                    help="regex stripped before the normalised join")
    ap.add_argument("--min-cluster-n", type=int, default=5)
    ap.add_argument("--residual-name", default="RESIDUAL")
    ap.add_argument("--residual-clonal-calls", action="store_true",
                    help="emit a keep-list for the residual group too (L2: normally off)")
    ap.add_argument("--out-assignments", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--out-keep-dir", required=True, type=Path)
    args = ap.parse_args()

    samples = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    if not samples:
        sys.exit(f"ERROR: {args.samples} is empty")
    labels = read_cluster_table(args.cluster_table)

    # Normalised lookup; first label wins if two raw IDs collapse to one key.
    norm_labels: dict[str, str] = {}
    for raw, cluster in labels.items():
        norm_labels.setdefault(normalise(raw, args.id_pattern), cluster)

    assignments: list[tuple[str, str]] = []
    n_exact = n_norm = 0
    for s in samples:
        if s in labels:
            assignments.append((s, labels[s]))
            n_exact += 1
        else:
            hit = norm_labels.get(normalise(s, args.id_pattern))
            if hit:
                assignments.append((s, hit))
                n_norm += 1
            else:
                assignments.append((s, args.residual_name))

    counts = Counter(c for _, c in assignments)
    n_residual = counts.get(args.residual_name, 0)

    print(f"[assign_ibd_clusters] cohort:            {len(samples)}")
    print(f"[assign_ibd_clusters] cluster table:     {len(labels)}")
    print(f"[assign_ibd_clusters] matched exact:     {n_exact}")
    print(f"[assign_ibd_clusters] matched normalised:{n_norm}"
          f"   (pattern={args.id_pattern!r})")
    print(f"[assign_ibd_clusters] unlabelled:        {n_residual}"
          f"  -> {args.residual_name}")

    args.out_keep_dir.mkdir(parents=True, exist_ok=True)
    args.out_assignments.parent.mkdir(parents=True, exist_ok=True)

    with args.out_assignments.open("w") as fh:
        fh.write("sample\tcluster\n")
        for s, c in assignments:
            fh.write(f"{s}\t{c}\n")

    eligible: list[str] = []
    with args.out_summary.open("w") as fh:
        fh.write("cluster\tn_samples\tibd_eligible\treason\n")
        for cluster, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            is_residual = cluster == args.residual_name
            if is_residual and not args.residual_clonal_calls:
                ok, reason = False, "residual group: not a homogeneous population (L2)"
            elif n < args.min_cluster_n:
                ok, reason = False, f"n < min_cluster_n ({args.min_cluster_n})"
            else:
                ok, reason = True, "ok"
            if ok:
                eligible.append(cluster)
            fh.write(f"{cluster}\t{n}\t{str(ok).lower()}\t{reason}\n")
            print(f"[assign_ibd_clusters]   {cluster:<14} n={n:<5} "
                  f"ibd={'yes' if ok else 'no':<3} {reason}")

    for cluster in eligible:
        keep = args.out_keep_dir / f"{cluster}.keep"
        with keep.open("w") as fh:
            for s, c in assignments:
                if c == cluster:
                    fh.write(f"{s}\t{s}\n")

    # Stale keep-lists from an earlier config would silently be picked up by a
    # glob downstream; drop any that are no longer eligible.
    for stale in args.out_keep_dir.glob("*.keep"):
        if stale.stem not in eligible:
            stale.unlink()
            print(f"[assign_ibd_clusters] removed stale keep-list {stale.name}")

    if not eligible:
        sys.exit("ERROR: no cluster is IBD-eligible; check --min-cluster-n and the "
                 "cluster table join")
    print(f"[assign_ibd_clusters] IBD-eligible clusters: {', '.join(eligible)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
