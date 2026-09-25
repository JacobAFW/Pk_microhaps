# 05_panel_select — candidate panel selection (Siegel NB2)

**Status:** rules and script **authored and linted 2026-09-25, never run.** The
include is live in `workflow/Snakefile` and the DAG resolves for all ten
algorithm × size combinations, but no panel has been produced or reviewed. The
exit criteria below are untested.

Deliberately **not** in `FINAL_TARGETS`: a bare `snakemake` must not start
building panels. Run it explicitly:

```bash
source envs/activate.sh
snakemake --cores 3 outputs/panel/v1_2026-09-25/panel_greedy_100.tsv   # the v1 deliverable
snakemake --cores 3 panel_select_all                                   # all 10 combinations
```

## What this stage does

Reproduces Siegel's `2_selection_from_microhaplotype_candidates.ipynb` for Pk.
Applies the tighter selection filters (3–10 SNPs, He ≥ 0.60), collapses
overlapping windows into independent loci, allocates markers per chromosome by
length, and selects panels of varying size by two algorithms.

## Inputs

- `outputs/scan/<scan_tag>/heterozygome.tsv` (Stage 04)
- `data/reference/strain_preDB_version/strain_A1_H.1.Icor.fasta.fai` — contig
  lengths for the length-proportional allocation. **Not** the PlasmoDB-67 `.fai`,
  which carries `LT*` accessions instead of the VCF's `ordered_PKNH_*_v2` naming
  (the same trap Stage 02b hit).

## Outputs

Version-tagged by `panel_select.version_tag`, so re-runs never overwrite:

- `outputs/panel/<tag>/candidates.<algorithm>_<size>.tsv` — selection-eligible loci
- `outputs/panel/<tag>/panel_<algorithm>_<size>.tsv` for `algorithm` in
  {greedy, evenly_spaced} and `size` in {50, 100, 150, 200, 250}

`outputs/panel/<tag>/panel_greedy_100.tsv` is the headline v1 candidate panel.

## Rules

| Rule | Purpose |
|---|---|
| `select_panel` | One panel, wildcards `{algorithm}` × `{size}`. |
| `panel_select_all` | Aggregate over all 10 combinations. Not a final target. |

---

## Design inputs from the Stage 04 run (2026-09-25)

These are measured facts from `runs/2026-09-25_stage04-window-scan/`, not
assumptions. Each carries the command that reproduces it.

### 1. Overlapping windows are not independent markers

Stage 04 windows are 200 bp wide at a 50 bp step, so neighbouring candidate
windows overlap by up to 75% and **describe the same piece of genome**.

| | count |
|---|---|
| selection-eligible windows (3–10 SNPs, He ≥ 0.60) | **19,674** |
| non-overlapping merged loci | **5,915** |
| windows per locus | 3.33 |
| headroom for a 100-marker panel | **59×**, not 197× |

```bash
python - <<'EOF'
import pandas as pd
w = pd.read_csv("outputs/scan/v1_2026-09-25/windows_all.tsv.gz", sep="\t")
e = w[(w.n_variants.between(3, 10)) & (w.het_raw >= 0.60)]
n = 0
for c, s in e.groupby("chrom"):
    end = -1
    for a, b in zip(s.sort_values("window_start").window_start,
                    s.sort_values("window_start").window_end):
        if a > end: n += 1; end = b
        else: end = max(end, b)
print(len(e), "windows ->", n, "merged loci")
EOF
```

**Consequence, already implemented:** `panel_select.py` merges to loci, keeps one
representative window per locus (highest `het_raw`, then most SNPs, then lowest
start — a total order, so deterministic), and enforces
`panel_select.min_spacing_bp` between selected markers. Both counts are printed.
`min_spacing_bp: 10000` is a **starting value, not a locked decision** — it has
never been run, and it is the first thing to sanity-check.

### 2. The panel is regionally tuned, and it is now measurable

Only **17.17%** of the 51,246 discovery SNPs clear MAF ≥ 0.10 in *all three*
labelled ADMIXTURE clusters (8,801 SNPs). Within the smallest cluster,
Peninsular (n=35), **38.7%** of discovery SNPs are monomorphic.

```bash
source envs/activate.sh
python runs/2026-09-25_stage04-window-scan/phase03_structure.py
# -> "SNPs common in ALL three clusters : 8,801 (17.17% of the pooled set)"
# -> within_cluster_maf.tsv, pct_monomorphic column
```

⚠ **Provenance:** the 17.17% (8,801 SNPs) is the **three-way intersection**, and
it exists only in `phase03_structure.py`'s stdout and in `phase_03_report.md`.
It is **not** recoverable from the TSVs — `within_cluster_maf.tsv` holds
per-cluster counts and `structure_attribution.tsv` holds the **union** (48,997),
never the intersection. Re-run the command above to regenerate it.

**Consequence, NOT implemented — needs a ruling.** A marker chosen on *pooled*
frequency may be near-uninformative in one cluster. Stage 05 should either
select on within-cluster frequency, or report per-cluster informativeness per
marker. Which of those depends on whether the v1 panel is Sabah/Mf-tuned
(defensible, and arguably right for the surveillance use case) or must work
across all three. **That is a NEEDS-JACOB.**

### 3. He is inflated by missingness, and switching statistic is cheap

`het_raw` counts haplotypes containing a missing call as distinct alleles
(decision L1, faithful to Siegel). At 6.5% per-call missingness this inflates He,
most in the high-SNP-count windows selection prefers. Recomputing the
heterozygome on `het_nomiss` instead costs **637 of 20,544 candidates (3.10%)**.

The Stage 04 recommendation is to keep `het_raw` as the Siegel-faithful filter
and carry `het_nomiss` as a **tiebreaker** when two loci are otherwise equal.
`panel_select.py` already writes `het_nomiss` into every panel row so the
tiebreaker is one line away. **Open decision before the panel is locked.**

## Known rough edges in `panel_select.py` — fix at the first commissioned run

The script has been smoke-tested at function level against the real heterozygome
(`runs/2026-09-25_stage04-window-scan/smoke_test_panel_select.py` — both
algorithms yield 100 markers, spacing respected, no all-NaN column), but it has
never been run end to end through Snakemake. Four things were deliberately left
alone rather than changed after review, because altering behaviour with no run to
validate it would trade inspected code for uninspected code:

1. **One real bug was found and fixed at review.** `pick_greedy` sorted on
   `n_variants`, which `main()` renames to `n_snps` before selection — so every
   greedy call, including the v1 deliverable, raised `KeyError`. `py_compile`,
   `snakemake --lint` and `snakemake -n` all passed on it. Fixed; the smoke test
   now covers it.
2. **`best_per_locus` uses `groupby(...).first()`**, which takes the first
   *non-null* value per column rather than the first row. Harmless today (the
   heterozygome has no NaNs, confirmed), but it would silently splice two rows
   together if Stage 04 ever emitted one. `.drop_duplicates("locus_id")` is the
   safer idiom.
3. **A shortfall exits 0.** If `min_spacing_bp` is binding, the script warns on
   stdout and still writes a file named `panel_<algo>_<N>.tsv` containing fewer
   than N markers, with no record of the shortfall inside the file. Prefer a
   non-zero exit, or a `requested_size` column, so the evidence is on disk rather
   than only in a log.
4. **`merge_to_loci`'s docstring says "overlapping/bookending"**, but the
   condition `if s > end` keeps exactly-touching windows separate. That is what
   produces the reported 5,915; merging bookends would give 5,630. The number is
   right, the wording is loose.

## Exit criterion (untested)

A 100-marker panel exists with an He distribution roughly comparable to Siegel's
Pv exemplar (median He 0.70–0.81 across regions). Inter-marker spacing is
reasonably uniform; no chromosome dropped entirely.

Given the Stage 04 heterozygome's `het_raw` median of 0.8058, a panel selected
greedily will sit **above** Siegel's band. That is expected — see the Stage 04
chapter on why the Pk discovery set is denser — and should be reported against
`het_nomiss` as well, so the comparison is not resting on the inflated statistic.

## Parameters (`workflow/config.yaml`, `panel_select:`)

| Key | Value | Meaning |
|---|---|---|
| `snp_min_per_window` | 3 | selection floor |
| `snp_max_per_window` | 10 | diminishing returns above 10; 2^10 = 1024 combos |
| `heterozygosity_min` | 0.60 | selection floor, on `het_raw` (L1) |
| `panel_sizes` | 50–250 | one panel per size |
| `algorithms` | greedy, evenly_spaced | |
| `default_panel_size` | 100 | the v1 deliverable |
| `min_spacing_bp` | 10000 | **new**, see design input 1. Starting value. |
| `version_tag` | v1_2026-09-25 | output subdirectory — bump, never overwrite |

## Blocking decisions before this stage is commissioned

1. **The 51k-SNP ruling** — whether MAF 0.10 stands. Evidence and a
   recommendation (keep it) are in the Stage 04 chapter and
   `runs/2026-09-25_stage04-window-scan/phase_03_report.md`. Changes the input.
2. **Cross-cluster vs Mf-tuned** (design input 2). Changes the selection criterion.
3. **`het_raw` vs `het_nomiss`** (design input 3). Changes the ranking.
4. **`min_spacing_bp`** (design input 1). Changes marker independence.

Stage 05 will run without these being settled, but the panel it produces would be
provisional in four separate ways.
