#!/usr/bin/env bash
# import_popgen_metadata.sh — copy the sample metadata table from Pop-gen
#
# Usage:
#   bash scripts/sh/import_popgen_metadata.sh <out_tsv>
#
# Resolution order:
#   1. POPGEN_PIPELINE_PATH env var, if set (preferred).
#   2. Sibling default: ../Pop-gen_pipeline relative to this project root.
#
# Source file priority (first match wins):
#   $POPGEN_PATH/data/metadata/fws_meta_indonesia_with_state.tsv
#   $POPGEN_PATH/outputs/setup/fws_meta_indonesia_with_state.tsv
#   $POPGEN_PATH/data/metadata/samples.tsv
#
# Output schema (TSV):
#   sample <TAB> country <TAB> region_geo
# where region_geo collapses Indonesia provinces + Malaysian states +
# any other country/state combos into the geographic group label used by
# Stage 04 onwards (Siegel's per-region heterozygosity reporting).
#
# CONTRACT
#   - The Pop-gen metadata file uses 'fws_meta_indonesia_with_state.tsv' by
#     convention. Its column names are not fully stable across re-runs; this
#     script auto-detects the relevant columns from the header.
#   - If the source file is missing, write a header-only stub. Stage 02's
#     merge_sample_qc.py will then abort with a helpful error.

set -euo pipefail

OUT=${1:-}
if [[ -z "$OUT" ]]; then
    echo "Usage: $0 <out_tsv>" >&2
    exit 1
fi
mkdir -p "$(dirname "$OUT")"

PROJECT_ROOT="$( cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd )"
POPGEN_DEFAULT="$( cd "$PROJECT_ROOT/.." 2>/dev/null && pwd )/Pop-gen_pipeline"
POPGEN_PATH="${POPGEN_PIPELINE_PATH:-$POPGEN_DEFAULT}"

CANDIDATES=(
    "$POPGEN_PATH/data/metadata/fws_meta_indonesia_with_state.tsv"
    "$POPGEN_PATH/outputs/setup/fws_meta_indonesia_with_state.tsv"
    "$POPGEN_PATH/data/metadata/samples.tsv"
)

SRC=""
for c in "${CANDIDATES[@]}"; do
    if [[ -s "$c" ]]; then
        SRC="$c"
        break
    fi
done

if [[ -z "$SRC" ]]; then
    cat <<EOF >&2

==> Pop-gen metadata table NOT found.
    Tried:
$(printf "      %s\n" "${CANDIDATES[@]}")

    Either:
      (a) re-run envs/install.sh after putting Pop-gen_pipeline next to this
          project, or set POPGEN_PIPELINE_PATH and re-run, or
      (b) drop a hand-curated TSV at data/metadata/samples.tsv with columns
          sample, country, region_geo.

    Writing a header-only stub; Stage 02's merge_sample_qc.py will abort with
    a clear error until this is resolved.

EOF
    printf "sample\tcountry\tregion_geo\n" > "$OUT"
    exit 0
fi

echo "==> Importing sample metadata from $SRC"

# Auto-detect columns. Sample-ID column lives under 'sample' or 'Sample' or
# 'sample_id'. Country in 'country' / 'Country'. State/region in any of:
# 'state', 'State', 'Region', 'region', 'province', 'Province'.
python3 - "$SRC" "$OUT" <<'PYEOF'
import csv
import re
import sys
from pathlib import Path

src = Path(sys.argv[1])
out = Path(sys.argv[2])

with src.open(newline="") as fh:
    sniffer = csv.Sniffer()
    sample = fh.read(4096)
    fh.seek(0)
    dialect = sniffer.sniff(sample, delimiters="\t,")
    reader = csv.DictReader(fh, dialect=dialect)
    headers = [h.strip() for h in reader.fieldnames or []]
    if not headers:
        sys.exit("Source file has no header row.")

    def find(*candidates):
        for h in headers:
            if h.lower() in {c.lower() for c in candidates}:
                return h
        return None

    sample_col  = find("sample", "sample_id", "ID", "Sample")
    country_col = find("country", "Country", "nation")
    region_col  = find("state", "State", "region", "Region",
                       "province", "Province", "region_geo")

    if not sample_col:
        sys.exit(f"Could not find a sample-ID column in {headers!r}")

    rows = list(reader)

with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "country", "region_geo"])
    for r in rows:
        sid = (r.get(sample_col) or "").strip()
        if not sid:
            continue
        country = (r.get(country_col) or "").strip() if country_col else ""
        region  = (r.get(region_col)  or "").strip() if region_col  else ""
        # region_geo := state/region if present, else country
        region_geo = region if region else country
        # Canonicalise: replace whitespace with underscore for downstream joins.
        region_geo = re.sub(r"\s+", "_", region_geo)
        w.writerow([sid, country, region_geo])

n = sum(1 for _ in open(out)) - 1
print(f"==> Wrote {n} samples to {out}", file=sys.stderr)
PYEOF
