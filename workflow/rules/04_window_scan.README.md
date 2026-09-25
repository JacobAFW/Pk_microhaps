# 04_window_scan — Pk heterozygome via sliding-window scan (Siegel NB1 phase 3)

**Status:** implemented 2026-09-25. Replaces the earlier stub, whose planned
`outputs/scan/all_windows.tsv` naming and 1,000–5,000-window exit criterion both
contradicted `handoffs/04_window_scan.md`. The handoff is the authority; the stub
text is superseded and is not a second opinion.

## What this stage does

Slides a 200 bp window in 50 bp steps along the Stage 03 discovery SNP set, one
chromosome at a time, and scores every window containing at least one SNP. The
scored table is filtered to the **heterozygome** — the candidate windows Stage 05
selects the panel from.

Ported from `external/clones/vivax-mhaps/notebooks/1_evaluating_marker_candidates.ipynb`
cells 36–41. The notebook's *markdown* says windows with ">1 variant" are stored;
its *code* keeps `n_variants != 0`. This follows the code, so single-SNP windows
are scored and written, and the ≥3-SNP requirement is applied only at the
heterozygome cut.

Per window:

| Statistic | Definition |
|---|---|
| `n_variants` | SNPs inside the window |
| `n_alleles` | distinct sample haplotype columns (`np.unique(gt, axis=1)`) |
| `het_raw` | `1 − Σf²` over those columns — Siegel's He |
| `entropy` | `−Σ f·ln f` |
| `n_alleles_with_missing` / `n_alleles_with_het` | Siegel's flag columns |
| `het_nomiss`, `entropy_nomiss`, `n_alleles_nomiss`, `n_samples_used` | the same statistics after haploidifying and dropping haplotypes with any missing call |
| `delta_het` | `het_raw − het_nomiss` |
| `frac_in_core` | fraction of the 200 bp span inside `pk_core_coding.bed` |

## Locked decisions

**L1 — missing-call policy: compute both.** `het_raw` counts missing-bearing
haplotypes as distinct alleles, which is what NB1 does (it only *flags* them, via
`unique_alleles_with_missing_index`). `het_nomiss` drops them before computing
frequencies. **The heterozygome is filtered on `het_raw` only**, so the scan stays
faithful and comparable to the paper, but `het_nomiss` and `delta_het` ride along
on every row of both outputs. This matters here in a way it did not for Pv: at
6.5% per-call missingness (deviation D1) the chance a window haplotype carries a
missing call is ~18% at 3 SNPs, ~29% at 5 and ~49% at 10 — so the inflation is
biased toward exactly the high-SNP-count windows panel selection prefers.

**L2 — het-call policy: match Siegel, and apply it to one branch only.**
`het_raw`, `entropy`, `n_alleles`, `n_alleles_with_missing` and
`n_alleles_with_het` are computed on the **un-haploidified diploid** sample
columns, which is what NB1 cell 36 does. `haploidify_samples()` with Siegel's
seed **250523** is applied **only** before the `het_nomiss` branch. Haploidifying
globally would make `het_raw` non-comparable to the paper and `delta_het`
meaningless. The seed is written into `window_scan_summary.tsv`. Impact is small
either way — the discovery VCF's het rate is 0.48% (deviation D2) — but the
policy is recorded rather than implicit.

One haploid realisation is drawn per chromosome after a single `np.random.seed`,
with chromosomes walked in sorted order, so the draw is reproducible and is
identical for every window overlapping a given variant.

**L3 — window placement vs core boundaries: report, don't filter (v1).** Siegel
anchors windows on core-restricted variants but never requires the 200 bp span to
sit inside a core interval. `pk_core_coding.bed` is 12,006 CDS fragments, so many
windows straddle a boundary — and an amplicon spanning a masked region is a real
primer-design problem later. `frac_in_core` is therefore reported on every row
and **never filtered on**; the cut is Stage 05's call.

**L4 — window de-duplication: replicate Siegel exactly.** Windows carrying an
identical set of variant positions collapse to the first occurrence. This is
NB1's `np.unique(positions, return_index=True)` block, it is easy to miss, and at
a 50 bp step over sparse SNPs it removes a large fraction of windows. `dedupe_windows()`
carries its own self-test, and the summary reports the count **both before and
after** the collapse.

## Inputs

| Path | From |
|---|---|
| `outputs/qc/snps.discovery.vcf.gz` (+ `.csi`) | Stage 03 |
| `outputs/qc/snps.summary.tsv` | Stage 03 — variant-count guard |
| `outputs/qc/discovery_samples.declonal.txt` | Stage 02b — cohort guard |
| `data/reference/pk_core_coding.bed` | Stage 01 — `frac_in_core` |

## Outputs

Version-tagged by `window_scan.version_tag` so re-runs never overwrite:

- `outputs/scan/<tag>/windows_all.tsv.gz` — every unique non-empty window, 17 columns
- `outputs/scan/<tag>/heterozygome.tsv` — the subset with `n_variants ≥ 3` **and** `het_raw ≥ 0.50`
- `outputs/scan/<tag>/window_scan_summary.tsv` — `scope / metric / value`: parameters, the funnel, per-chromosome counts, and distribution medians
- `reports/figures/<figures_subdir>/{window_het_distribution, window_het_vs_entropy, window_manhattan_het, window_het_raw_vs_nomiss}.{pdf,png}`

Column order in both tables is fixed by the handoff contract:
`chrom, window_start, window_end, midpoint, n_variants, variant_positions,
n_alleles, het_raw, entropy, n_alleles_with_missing, n_alleles_with_het,
n_alleles_nomiss, het_nomiss, entropy_nomiss, delta_het, n_samples_used,
frac_in_core`.

## Rules

| Rule | Purpose |
|---|---|
| `window_scan` | The scan. Writes all three tables. Single-threaded. |
| `plot_window_scan` | Four figures, each as vector PDF + 300 dpi PNG. |

## Parameters (`workflow/config.yaml`, `window_scan:`)

| Key | Value | Meaning |
|---|---|---|
| `size_bp` | 200 | max Illumina amplicon size |
| `step_bp` | 50 | 75% overlap |
| `min_snps` | 3 | heterozygome floor |
| `min_heterozygosity` | 0.50 | heterozygome floor, on `het_raw` (L1) |
| `seed` | 250523 | Siegel's haploidify seed (L2) |
| `missing_policy` | both | L1; recorded in the summary |
| `version_tag` | `v1_2026-09-25` | output subdirectory — bump, never overwrite |
| `figures_subdir` | `04_window_scan` | figure subdirectory under `reports/figures/` |

## Guards

Four, all fail loud:

1. **No overwrite.** The script exits if any of its three outputs already exists.
2. **Cohort.** VCF sample set must equal `discovery_samples.declonal.txt` — the
   anti-silent-loss check at the Stage 03 → Stage 04 boundary (expected 593).
3. **Variant count.** Variants read must equal `snps.summary.tsv`, overall *and*
   per chromosome (expected 51,251).
4. **L4 self-test.** `dedupe_windows()` asserts keep-first on a fixture before
   any real data is touched.

## Exit criteria

Siegel's Pv funnel is the scaling anchor: 13,084 SNPs → 13,498 unique windows →
3,830 at He ≥ 0.5 → 1,110 at 3–10 SNPs & He ≥ 0.6. Our 51,251 SNPs is ~3.9×
that, so *proportional* scaling predicts roughly 45–55K unique windows, ~12–16K
at He ≥ 0.5 and ~3.5–4.5K final candidates. **That is an expectation to test, not
a target to hit** — the point of this stage is to find out whether the SNP excess
(deviation C2) is real Pk diversity or a leaky core-genome proxy (B1).

Pass:

- `heterozygome.tsv` non-empty, funnel written to the summary
- the He-by-variant-count boxplot reproduces Siegel's shape: monotone rise,
  diminishing returns past ~10 variants
- no chromosome holds a candidate-window share wildly out of line with its SNP
  share in `snps.summary.tsv`
- `delta_het` reported and its magnitude stated in the chapter

## Stop conditions

Stop and escalate if any of these fires:

- top-He windows dominated by *SICAvar* / *kir* → the core BED is leaking
- heterozygome empty, or fewer than ~500 windows pass
- unique windows above ~80K or below ~20K → check L4 first
- `delta_het` median above ~0.05 → missingness is distorting selection more than
  L1 assumed, and L1 should be revisited before Stage 05

## Known gaps carried forward

- **Per-region He is not reported.** 912/990 samples carry `region_geo='UNKNOWN'`
  (deviation A6); the region-metadata strategy is an open decision. The chapter
  is built so a per-region breakdown drops in later, and degrades to a single
  `UNKNOWN` stratum until then. Region is never imputed.
- **The core-BED hypervariable guarantee is not yet a pipeline rule.**
  `scripts/py/audit_core_coding_bed.py` exists and `DECISIONS.md` describes it as
  an enforced Stage 01 rule, but no such rule is wired. It is run standalone for
  this stage's sanity check; wiring it belongs to the next run.

## Reference implementation

`external/clones/vivax-mhaps/notebooks/1_evaluating_marker_candidates.ipynb`,
cells 36–41. Read the code, not the markdown.
