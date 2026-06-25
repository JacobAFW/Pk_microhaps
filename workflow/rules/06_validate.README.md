# 06_validate — panel evaluation (Siegel NB2 + paper Methods)

**Status:** stub — rules not yet written.

## What this stage will do
Replicate Siegel's panel-evaluation steps for the candidate Pk 100-marker panel:
- `paneljudge` IBD CI/RMSE simulations across r ∈ {0, 0.25, 0.5, 0.75, 1.0}
- per-region HE + effective-cardinality table
- (deferred to v1.5) hmmIBD on real Pop-gen sample pairs
- (deferred to v1.5) PCoA on microhap-based pairwise distance
- (deferred to v1.5) BALK-classifier-style geographic origin prediction

v1 ships with paneljudge simulation only — that's the headline scientific claim ("RMSE < 0.1 at r=0.5 with 100 microhaps") and is fast, deterministic, and reference-only. Real-data validation is moved to v1.5.

## Planned inputs
- `outputs/panel/panel_greedy_100.tsv` (from Stage 05)
- `outputs/qc/snps.discovery.vcf.gz` (from Stage 03)
- `outputs/qc/sample_metadata.tsv` (from Stage 02; needs region column)

## Planned outputs
- `outputs/validation/paneljudge_per_region.tsv` — per-region CI width and RMSE at each r
- `outputs/validation/panel_he_per_region.tsv` — He + effective cardinality per region
- `reports/figures/paneljudge_ci_per_region.png` — Fig. 3 equivalent
- `reports/figures/paneljudge_rmse_vs_panel_size.png` — Supp. Fig. 2 equivalent (across panel sizes)

## Planned scripts
- `scripts/R/paneljudge_eval.R` — argv-driven R script wrapping `paneljudge::simulate_ibd()` + RMSE/CI computation
- `scripts/R/panel_diversity_summary.R` — per-region He + effective cardinality

## Exit criterion (v1 acceptance)
- 100-marker panel achieves RMSE < 0.1 at r=0.5 in the Pk regional groups with the largest sample counts (Sabah, Indonesia provinces). For groups with small N, expect wider CIs — flag rather than treating as failure.
- Per-region median HE comparable to Siegel's Pv ranges (≥ 0.6 floor in all groups; aspirational 0.7+).

## Tools (installed by `INSTALL_VALIDATION=1 bash envs/install.sh`)
- `paneljudge` (R, github.com/aimeertaylor/paneljudge)
- `hmmIBD` v3 — built from source; deferred to v1.5
- `THE REAL McCOIL` v2 — built from source; deferred to v1.5
