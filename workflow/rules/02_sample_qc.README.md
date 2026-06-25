# 02_sample_qc — monoclonal sample filter (Siegel NB1 phase 1)

**Status:** implemented (HANDOFF #03, 2026-05-12). Not yet run.

## What this stage does
Filters the 990-sample Pk cohort to high-quality monoclonal samples for microhap discovery, using Siegel et al's criteria:
- **Fws ≥ 0.95** (monoclonal proxy) — sourced from Pop-gen Stage 2 (`outputs/setup/fws_from_popgen.tsv`).
- **callable_fraction ≥ 0.50** — fraction of SNP sites within `pk_core_coding.bed` with a non-missing genotype. This is a v1 proxy for Siegel's coverage-based callable; conservative when paired with Fws ≥ 0.95.

Sample metadata (country, state/region) is imported from Pop-gen's `fws_meta_indonesia_with_state.tsv` so downstream stages have a `region_geo` label (country/state per the decision locked in MEMORY.md).

## Inputs
- `data/vcf/merged_popgen.vcf.gz` (+ .csi) — Stage 01
- `data/reference/pk_core_coding.bed` — Stage 01
- `outputs/setup/fws_from_popgen.tsv` — Stage 01
- `data/metadata/samples.tsv` — written by `import_popgen_metadata` (this stage)

## Outputs
- `data/metadata/samples.tsv` — sample, country, region_geo.
- `outputs/qc/callable_fraction.tsv` — sample, n_called, n_total, callable_fraction.
- `outputs/qc/sample_metadata.tsv` — joined per-sample table with the include flag.
- `outputs/qc/discovery_samples.txt` — one sample ID per line; input to Stage 03.
- `reports/figures/fws_distribution.png`
- `reports/figures/callable_fraction_distribution.png`

## Rules
| Rule | Purpose |
|---|---|
| `import_popgen_metadata` | Pulls country/state metadata from Pop-gen. |
| `compute_callable_fraction` | Per-sample non-missing-GT fraction over `pk_core_coding.bed`. bcftools query streaming. |
| `merge_sample_qc` | Joins Fws + callable + metadata; applies thresholds; writes the discovery list. Hard-asserts on Fws ∈ [0,1] and ≥50 samples passing. |
| `plot_sample_qc` | Sanity plots, coloured by pass/fail. |

## Exit criterion
≥50 samples pass both thresholds. Plot looks bimodal-ish (clear monoclonal/polyclonal separation). Siegel's Pv ratio was 615 / 1,816 ≈ 34%; we'd hope for similar or better against the 990-sample Pk input.

## Open caveats / TODOs
- **Fws/Proportion equivalence.** The Pop-gen value column is literally `Proportion` (moimix convention). Stage 02 assumes this equals Fws and asserts it's in [0,1]. Worth a 30-second confirmation with the Pop-gen author. Logged in CONTEXT-software.md → Open issues #2.
- **Coverage-based callable.** v1 proxy is SNP-call-rate over core+coding. If BAMs become accessible, swap in true per-base coverage callable.
- **Region labels.** Decided 2026-05-12: use country/state from metadata (region_geo). ADMIXTURE-cluster alternative remains parked in PLAN.md if Stage 04 results suggest the geographic grouping isn't carving the cohort well.
