#!/usr/bin/env python
"""
annotate_candidate_windows.py — which genes do the heterozygome's candidate
windows sit in?

`handoffs/04_window_scan.md` requires the top 20 genes by candidate-window count
("it is cheap and it is the kind of table a reviewer will ask for"), and the
narrowed spatial-clustering diagnostic needs the gene-family composition of the
top-He windows: the *SICAvar* / *kir* question is closed (deviation B1), so if
the excess windows are spatially clustered the remaining job is to name which
genes they are clustered in.

Read-only on the GFF and read-only on the heterozygome. **No filter is applied to
the heterozygome** — every candidate window is annotated, and windows overlapping
no annotated gene are counted and reported rather than dropped, so the per-window
table reconciles row-for-row with `heterozygome.tsv`.

A window is assigned to every gene whose span it overlaps, so one window can
appear against more than one gene (genes overlap and nest in this annotation).
The per-gene window counts therefore sum to at least the number of annotated
windows, never to the heterozygome row count — both totals are printed.

Gene labels come from GFF attributes only. `Name=` when present, otherwise the
bare `ID=`; the `labelled_by` column says which, per gene, so an unnamed locus is
visible as such rather than silently blank.

Flagged families use the same classifier as `scripts/py/audit_core_coding_bed.py`
so the two agree by construction: *SICAvar*, *kir*, and pseudogenes.

Usage:
    python scripts/py/annotate_candidate_windows.py \
        --heterozygome outputs/scan/v1_2026-09-25/heterozygome.tsv \
        --gff data/reference/pk_annotation.lifted.gff \
        --top-n 20 \
        --out-top outputs/scan/v1_2026-09-25/candidate_genes_top20.tsv \
        --out-windows outputs/scan/v1_2026-09-25/candidate_window_genes.tsv
"""

import argparse
import collections
import re
import sys
from pathlib import Path

import pandas as pd

# Gene-level feature types in the lifted GFF. Exons/CDS/mRNA are children and
# would double-count, so they are deliberately excluded.
GENE_TYPES = ("protein_coding_gene", "pseudogene", "ncRNA_gene", "gene")


def classify(name, description, feature_type):
    """Flagged family for a gene, or 'unflagged'. Same rules as
    scripts/py/audit_core_coding_bed.py, so the two reconcile."""
    text = f"{name} {description}".lower()
    if "sicavar" in text or re.search(r"\bsica\b", text):
        return "SICAvar"
    if re.search(r"\bkir\b", text):
        return "kir"
    if feature_type == "pseudogene":
        return "pseudogene"
    return "unflagged"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--heterozygome", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--out-top", required=True)
    p.add_argument("--out-windows", required=True)
    return p.parse_args()


def load_genes(gff_path):
    """chrom -> list of (start, end, gene_id, label, labelled_by, family, biotype).
    Coordinates stay 1-based inclusive, as in the GFF and as in the window table."""
    genes = collections.defaultdict(list)
    with Path(gff_path).open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            attrs = f[8]
            gid = (re.search(r"ID=([^;]*)", attrs) or [None, ""])[1]
            name = (re.search(r"Name=([^;]*)", attrs) or [None, ""])[1]
            desc = (re.search(r"description=([^;]*)", attrs) or [None, ""])[1]
            biotype = (re.search(r"ebi_biotype=([^;]*)", attrs) or [None, ""])[1]
            label, how = (name, "Name") if name else (gid, "ID")
            genes[f[0]].append((int(f[3]), int(f[4]), gid, label, how,
                                classify(name, desc, f[2]), biotype or f[2]))
    for c in genes:
        genes[c].sort()
    return genes


def overlapping_genes(genes_on_chrom, start, end):
    """Every gene whose 1-based inclusive span intersects [start, end].

    Linear scan per window: genes nest and overlap here, so an interval-tree
    shortcut would need care for no real gain — the heterozygome is thousands of
    rows against ~5,500 genes per genome.
    """
    return [g for g in genes_on_chrom if g[0] <= end and g[1] >= start]


def main():
    args = parse_args()

    het = pd.read_csv(args.heterozygome, sep="\t")
    if het.empty:
        sys.exit(f"ERROR: {args.heterozygome} has no rows — nothing to annotate.")
    for col in ("chrom", "window_start", "window_end", "het_raw", "n_variants"):
        if col not in het.columns:
            sys.exit(f"ERROR: {args.heterozygome} lacks required column {col!r}.")

    genes = load_genes(args.gff)
    n_genes = sum(len(v) for v in genes.values())
    unnamed = sum(1 for v in genes.values() for g in v if g[4] == "ID")
    print(f"==> {n_genes:,} gene-level features from {args.gff} "
          f"({len(genes)} contigs; {unnamed:,} carry no Name= and are labelled by ID)")
    print(f"==> {len(het):,} candidate windows from {args.heterozygome} "
          f"(no filter applied)")

    per_window, per_gene = [], collections.defaultdict(
        lambda: {"n_windows": 0, "max_het_raw": float("-inf")})
    n_no_gene = 0
    for r in het.itertuples():
        hits = overlapping_genes(genes.get(r.chrom, []), r.window_start, r.window_end)
        if not hits:
            n_no_gene += 1
        per_window.append({
            "chrom": r.chrom,
            "window_start": r.window_start,
            "window_end": r.window_end,
            "n_variants": r.n_variants,
            "het_raw": r.het_raw,
            "n_genes": len(hits),
            "gene_ids": ",".join(g[2] for g in hits),
            "gene_labels": ",".join(g[3] for g in hits),
            "families": ",".join(sorted({g[5] for g in hits})),
        })
        for g in hits:
            key = (r.chrom, g[2], g[3], g[4], g[5], g[6], g[0], g[1])
            per_gene[key]["n_windows"] += 1
            per_gene[key]["max_het_raw"] = max(per_gene[key]["max_het_raw"], r.het_raw)

    pw = pd.DataFrame(per_window)
    Path(args.out_windows).parent.mkdir(parents=True, exist_ok=True)
    pw.to_csv(args.out_windows, sep="\t", index=False, float_format="%.6g")

    rows = [{
        "chrom": k[0], "gene_id": k[1], "gene_label": k[2], "labelled_by": k[3],
        "family": k[4], "biotype": k[5], "gene_start": k[6], "gene_end": k[7],
        "gene_length_bp": k[7] - k[6] + 1,
        "n_candidate_windows": v["n_windows"], "max_het_raw": v["max_het_raw"],
    } for k, v in per_gene.items()]
    gene_df = (pd.DataFrame(rows)
               .sort_values(["n_candidate_windows", "max_het_raw"], ascending=False)
               .reset_index(drop=True))
    gene_df.head(args.top_n).to_csv(args.out_top, sep="\t", index=False,
                                    float_format="%.6g")

    # Reconciliation, printed so the report can cite it without recomputing.
    annotated = len(het) - n_no_gene
    assignments = int(gene_df["n_candidate_windows"].sum())
    print(f"\n==> candidate windows          : {len(het):,}")
    print(f"    overlapping >=1 gene       : {annotated:,}")
    print(f"    overlapping NO gene        : {n_no_gene:,} "
          f"({100 * n_no_gene / len(het):.2f}%)")
    print(f"    genes with >=1 window      : {len(gene_df):,}")
    print(f"    window-gene assignments    : {assignments:,} "
          f"(>= annotated windows, because windows can span several genes)")
    fam = collections.Counter()
    for k, v in per_gene.items():
        fam[k[4]] += v["n_windows"]
    print(f"    assignments by family      : "
          + ", ".join(f"{k} {v:,}" for k, v in fam.most_common()))
    print(f"\n==> wrote {args.out_top}\n            {args.out_windows}")


if __name__ == "__main__":
    main()
