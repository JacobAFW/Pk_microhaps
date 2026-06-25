#!/usr/bin/env bash
# build_core_coding_bed.sh — derive Pk core+coding BED for the window scan
#
# Usage:
#   bash scripts/sh/build_core_coding_bed.sh <gff> <mask_list> <out_bed>
#
# Strategy:
#   1. Pull CDS features from the GFF → BED (chrom, start-1, end).
#   2. Subtract the regions in mask_list (chrom:start-end format).
#   3. Apply chromosome-name translation if data/reference/chrom_translation.tsv
#      exists (CDS GFF uses ENA LT* accessions; VCF + mask use ordered_PKNH_*_v2
#      naming — the two need reconciling before bedtools can intersect cleanly).
#
# KNOWN ISSUE — see CONTEXT-software.md "Open issues" section. The PlasmoDB-67/68
# GFF chromosomes are named differently from the VCF/mask. This script will
# WARN loudly if a chr-translation table is missing and produce a BED with
# whichever naming the GFF uses. Stage 03 will fail on that BED until names
# are reconciled. Decision pinned in PLAN.md → Open decisions.

set -euo pipefail

GFF=${1:-}
MASK_LIST=${2:-}
OUT_BED=${3:-}

if [[ -z "$GFF" || -z "$MASK_LIST" || -z "$OUT_BED" ]]; then
    echo "Usage: $0 <gff> <mask_list> <out_bed>" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUT_BED")
mkdir -p "$OUT_DIR"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ---------------------------------------------------------------------------
# 1) GFF → CDS BED
#    GFF columns: seqid source feature start end score strand frame attrs
#    Filter to feature == CDS, convert to 0-based BED.
# ---------------------------------------------------------------------------
echo "==> Extracting CDS features from $GFF"
awk 'BEGIN{OFS="\t"} !/^#/ && $3 == "CDS" { print $1, $4 - 1, $5 }' "$GFF" \
  | sort -k1,1 -k2,2n \
  > "$TMP/cds.raw.bed"
n_cds=$(wc -l < "$TMP/cds.raw.bed")
echo "    $n_cds CDS intervals"

# Merge overlapping CDS intervals (alternative isoforms etc.)
if command -v bedtools &>/dev/null; then
    bedtools merge -i "$TMP/cds.raw.bed" > "$TMP/cds.merged.bed"
else
    # Lightweight awk merge if bedtools isn't available yet (pre-install)
    awk 'BEGIN{OFS="\t"; pc=""; ps=-1; pe=-1}
         {
           if ($1 == pc && $2 <= pe) {
             if ($3 > pe) pe = $3
           } else {
             if (pc != "") print pc, ps, pe
             pc = $1; ps = $2; pe = $3
           }
         }
         END { if (pc != "") print pc, ps, pe }' "$TMP/cds.raw.bed" \
         > "$TMP/cds.merged.bed"
fi
n_merged=$(wc -l < "$TMP/cds.merged.bed")
echo "    $n_merged merged CDS intervals"

# ---------------------------------------------------------------------------
# 2) Mask list (chrom:start-end) → BED
# ---------------------------------------------------------------------------
echo "==> Converting mask list to BED"
awk -F'[:-]' 'BEGIN{OFS="\t"} NF==3 { print $1, $2 - 1, $3 }' "$MASK_LIST" \
  | sort -k1,1 -k2,2n > "$TMP/mask.bed"
n_mask=$(wc -l < "$TMP/mask.bed")
echo "    $n_mask mask intervals"

# ---------------------------------------------------------------------------
# 3) Chromosome-name reconciliation (CRITICAL — see header)
# ---------------------------------------------------------------------------
TRANSLATION="data/reference/chrom_translation.tsv"
if [[ -f "$TRANSLATION" ]]; then
    echo "==> Applying chrom-name translation from $TRANSLATION"
    # Translate the CDS BED's chrom column (col 1) using the table.
    # Translation table format (TSV, header optional):
    #   gff_name<TAB>vcf_name
    awk -v tr="$TRANSLATION" '
        BEGIN {
            FS=OFS="\t"
            while ((getline line < tr) > 0) {
                if (line ~ /^#/ || line ~ /^gff_name/) continue
                split(line, a, "\t")
                if (length(a) >= 2) map[a[1]] = a[2]
            }
        }
        { if ($1 in map) $1 = map[$1]; print }' \
        "$TMP/cds.merged.bed" | sort -k1,1 -k2,2n > "$TMP/cds.translated.bed"
else
    echo "==> WARNING: $TRANSLATION not found." >&2
    echo "    GFF uses ENA accessions (e.g. LT727662); VCF + mask use" >&2
    echo "    ordered_PKNH_*_v2 naming. Stage 03 will fail on the resulting BED" >&2
    echo "    until a translation table is provided. See CONTEXT-software.md." >&2
    cp "$TMP/cds.merged.bed" "$TMP/cds.translated.bed"
fi

# ---------------------------------------------------------------------------
# 4) Subtract mask from CDS → core+coding BED
# ---------------------------------------------------------------------------
if command -v bedtools &>/dev/null; then
    bedtools subtract -a "$TMP/cds.translated.bed" -b "$TMP/mask.bed" > "$OUT_BED"
else
    echo "==> WARNING: bedtools not on PATH — writing un-subtracted BED." >&2
    echo "    Re-run after running envs/install.sh." >&2
    cp "$TMP/cds.translated.bed" "$OUT_BED"
fi

n_out=$(wc -l < "$OUT_BED")
total_bp=$(awk '{ s += $3 - $2 } END { print s }' "$OUT_BED")
echo "==> Wrote $OUT_BED  ($n_out intervals, $total_bp bp)"
