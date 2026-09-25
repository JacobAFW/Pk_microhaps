#!/usr/bin/env python
"""
audit_core_coding_bed.py — is `pk_core_coding.bed` actually excluding the
hypervariable multigene families?

Pv4 ships a curated core-genome annotation. Pk has none, so `pk_core_coding.bed`
is built as CDS minus Pop-gen's mapping-quality mask list. That mask was built for
pop-gen variant calling, not as a hypervariable-region annotation, so this script
provides the orthogonal check: locate annotated SICAvar / kir loci from the lifted
GFF and report how much of them survives into the BED and into the discovery SNP set.

Also emits a cleaned BED with any surviving flagged intervals removed, so the
result is enforced rather than merely reported.

Run standalone or from Snakemake (Stage 01).

Usage:
    python scripts/py/audit_core_coding_bed.py \
        --gff data/reference/pk_annotation.lifted.gff \
        --bed data/reference/pk_core_coding.bed \
        --vcf outputs/qc/snps.discovery.vcf.gz \
        --out-report outputs/setup/core_coding_audit.tsv \
        --out-bed data/reference/pk_core_coding.clean.bed \
        --max-flagged-snp-frac 0.01

Exit codes:
    0  audit passed (flagged-family SNP fraction within tolerance)
    2  audit FAILED — the mask is letting hypervariable families through.
       Do not proceed to panel selection; the core-genome proxy needs revisiting.

Baseline (2026-08-17, 990-sample merged_popgen VCF, 602-sample discovery cohort):
    raw CDS      10,279,998 bp   SICAvar 476,739 (4.64%)  kir 97,977 (0.95%)
    core BED      9,581,281 bp   SICAvar   6,312 (0.07%)  kir      0 (0.00%)
    discovery SNPs   49,803      SICAvar       0           kir      0
    -> 81% of everything the mask removed was SICAvar/kir CDS.
Any material departure from this on a new VCF or a new mask is a real signal.
"""

import argparse
import bisect
import collections
import gzip
import re
import sys
from pathlib import Path

# Gene families to flag. Pk's hypervariable multigene families — the analogue of
# the subtelomeric/hypervariable regions Pv4's curated core genome excludes.
GENE_TYPES = ("protein_coding_gene", "pseudogene", "ncRNA_gene", "gene")


def classify(name, description, feature_type):
    """Return the flagged family for a gene, or None."""
    text = f"{name} {description}".lower()
    if "sicavar" in text or re.search(r"\bsica\b", text):
        return "SICAvar"
    if re.search(r"\bkir\b", text):
        return "kir"
    if feature_type == "pseudogene":
        return "pseudogene"
    return None


def load_flagged_genes(gff_path):
    """chrom -> sorted [(start0, end, family)] for every flagged gene."""
    flagged = collections.defaultdict(list)
    with Path(gff_path).open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            attrs = f[8]
            name = (re.search(r"Name=([^;]*)", attrs) or [None, ""])[1]
            desc = (re.search(r"description=([^;]*)", attrs) or [None, ""])[1]
            fam = classify(name, desc, f[2])
            if fam:
                flagged[f[0]].append((int(f[3]) - 1, int(f[4]), fam))
    for chrom in flagged:
        flagged[chrom].sort()
    return flagged


class FamilyIndex:
    """Position -> flagged family lookup. Genes nest and overlap, so the
    backscan window must be generous; SELF-TEST below proves it is."""

    BACKSCAN = 30

    def __init__(self, flagged):
        self.idx = {c: ([s for s, _, _ in ivs], ivs) for c, ivs in flagged.items()}

    def at(self, chrom, pos0):
        if chrom not in self.idx:
            return None
        starts, ivs = self.idx[chrom]
        i = bisect.bisect_right(starts, pos0) - 1
        lo = max(0, i - self.BACKSCAN)
        for j in range(lo, min(len(ivs), i + 2)):
            s, e, fam = ivs[j]
            if s <= pos0 < e:
                return fam
        return None

    def self_test(self, flagged):
        """Every flagged gene's own midpoint must resolve to its own family.
        Guards against a silently-too-small backscan producing a clean-looking
        but meaningless audit."""
        total = misses = 0
        for chrom, ivs in flagged.items():
            for s, e, fam in ivs:
                total += 1
                if self.at(chrom, (s + e) // 2) != fam:
                    misses += 1
        return total, misses


def bed_composition(bed_path, index):
    """Per-family interval count and bp in a BED."""
    bp = collections.Counter()
    n = collections.Counter()
    rows = []
    with Path(bed_path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split()[:3]
            s, e = int(s), int(e)
            fam = index.at(c, s) or index.at(c, e - 1) or "unflagged"
            bp[fam] += e - s
            n[fam] += 1
            rows.append((c, s, e, fam))
    return bp, n, rows


def snp_composition(vcf_path, index):
    """Per-family SNP counts, overall and per chromosome."""
    total = collections.Counter()
    per_chrom = collections.defaultdict(collections.Counter)
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _ = line.split("\t", 2)
            fam = index.at(chrom, int(pos) - 1) or "unflagged"
            total[fam] += 1
            per_chrom[chrom][fam] += 1
    return total, per_chrom


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True, help="lifted GFF with gene names")
    p.add_argument("--bed", required=True, help="pk_core_coding.bed to audit")
    p.add_argument("--vcf", help="discovery SNP VCF (optional; skipped if absent)")
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-bed", help="write a BED with flagged intervals removed")
    p.add_argument("--max-flagged-snp-frac", type=float, default=0.01,
                   help="fail if flagged families hold more than this "
                        "fraction of discovery SNPs (default 0.01)")
    return p.parse_args()


def main():
    args = parse_args()

    flagged = load_flagged_genes(args.gff)
    index = FamilyIndex(flagged)

    total, misses = index.self_test(flagged)
    print(f"==> SELF-TEST: {total} flagged genes, {misses} midpoints misclassified")
    if misses:
        sys.exit(f"ERROR: position lookup is unreliable ({misses}/{total} misses). "
                 "Increase FamilyIndex.BACKSCAN. Audit results would be meaningless.")

    span = collections.Counter()
    for ivs in flagged.values():
        for s, e, fam in ivs:
            span[fam] += e - s
    print(f"    flagged gene span: "
          + ", ".join(f"{k} {v:,} bp" for k, v in sorted(span.items())))

    bed_bp, bed_n, rows = bed_composition(args.bed, index)
    bed_total = sum(bed_bp.values())
    print(f"\n==> {args.bed}: {sum(bed_n.values()):,} intervals / {bed_total:,} bp")
    for fam, v in bed_bp.most_common():
        print(f"    {fam:12s} {bed_n[fam]:6,} iv  {v:10,} bp  ({100*v/bed_total:5.2f}%)")

    report = [("metric", "family", "value", "pct")]
    for fam, v in bed_bp.items():
        report.append(("bed_bp", fam, str(v), f"{100*v/bed_total:.4f}"))
        report.append(("bed_intervals", fam, str(bed_n[fam]), ""))

    exit_code = 0
    if args.vcf and Path(args.vcf).exists():
        snp_total, per_chrom = snp_composition(args.vcf, index)
        n_snp = sum(snp_total.values())
        flagged_snp = sum(v for k, v in snp_total.items() if k != "unflagged")
        frac = flagged_snp / n_snp if n_snp else 0.0
        print(f"\n==> {args.vcf}: {n_snp:,} SNPs")
        for fam, v in snp_total.most_common():
            print(f"    {fam:12s} {v:7,}  ({100*v/n_snp:5.2f}%)")
        print(f"    flagged-family SNP fraction: {frac:.4%} "
              f"(tolerance {args.max_flagged_snp_frac:.2%})")
        for fam, v in snp_total.items():
            report.append(("snp_count", fam, str(v), f"{100*v/n_snp:.4f}"))
        for chrom in sorted(per_chrom):
            d = per_chrom[chrom]
            bad = sum(v for k, v in d.items() if k != "unflagged")
            report.append((f"snp_flagged:{chrom}", "any", str(bad),
                           f"{100*bad/sum(d.values()):.4f}"))
        if frac > args.max_flagged_snp_frac:
            print("\n*** AUDIT FAILED ***", file=sys.stderr)
            print(f"{frac:.2%} of discovery SNPs sit in SICAvar/kir/pseudogenes. "
                  "The mask is not excluding Pk hypervariable families. Panel "
                  "windows selected from these regions will look diverse and "
                  "will not amplify reliably. Stop before Stage 05.", file=sys.stderr)
            exit_code = 2
    else:
        print("\n==> No VCF given — BED-only audit.")

    with Path(args.out_report).open("w") as fh:
        for row in report:
            fh.write("\t".join(row) + "\n")
    print(f"\n==> Wrote {args.out_report}")

    if args.out_bed:
        kept = [r for r in rows if r[3] == "unflagged"]
        dropped = len(rows) - len(kept)
        with Path(args.out_bed).open("w") as fh:
            for c, s, e, _ in kept:
                fh.write(f"{c}\t{s}\t{e}\n")
        print(f"==> Wrote {args.out_bed} "
              f"({len(kept):,} intervals; dropped {dropped:,} flagged)")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
