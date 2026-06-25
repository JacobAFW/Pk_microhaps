#!/usr/bin/env bash
# import_popgen_fws.sh — copy the Fws table from Pop-gen Stage 2 into this project
#
# Usage:
#   bash scripts/sh/import_popgen_fws.sh <out_tsv>
#
# Resolution order:
#   1. POPGEN_PIPELINE_PATH env var, if set.
#   2. Sibling default: ../Pop-gen_pipeline relative to this project root.
#
# Source file: outputs/moi/fws_MOI.tsv inside Pop-gen.
# Schema: sample<TAB>Proportion (where "Proportion" is the moimix Fws-like statistic).
# Output: standardised TSV with columns sample<TAB>fws.
#
# If Pop-gen isn't reachable, emit a header-only stub and a clear stderr
# message. Stage 02 will then either (a) accept the stub and fall back to
# moimix, or (b) abort with a helpful error.

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
SRC="$POPGEN_PATH/outputs/moi/fws_MOI.tsv"

if [[ -f "$SRC" ]]; then
    echo "==> Importing Fws from $SRC"
    # Standardise: the Pop-gen file uses 'Proportion' as the value column;
    # rename to 'fws' for clarity here.
    awk 'BEGIN{OFS="\t"}
         NR == 1 { print "sample", "fws"; next }
         { print $1, $2 }' "$SRC" > "$OUT"
    n=$(($(wc -l < "$OUT") - 1))
    echo "==> Wrote $OUT ($n samples)"
else
    cat <<EOF >&2

==> Pop-gen Fws table NOT found.
    Tried: $SRC

    Either:
      (a) re-run envs/install.sh after putting Pop-gen_pipeline next to this
          project, or set POPGEN_PIPELINE_PATH and re-run, or
      (b) tell Stage 02 to recompute Fws via moimix (slower; needs a stricter
          MAF-filtered VCF as input).

    Writing a header-only stub so the rule still completes; Stage 02 will
    handle the empty case explicitly.

EOF
    printf "sample\tfws\n" > "$OUT"
fi
