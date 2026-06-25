# 05_panel_select — candidate panel selection (Siegel NB2)

**Status:** stub — rules not yet written.

## What this stage will do
Reproduce Siegel's `2_selection_from_microhaplotype_candidates.ipynb` for Pk. Apply the tighter selection filters (3–10 SNPs, He ≥ 0.6), allocate markers per chromosome by length, and select panels of varying size using two algorithms.

## Planned inputs
- `outputs/scan/heterozygome.tsv` (from Stage 04)
- `data/reference/pk_core_coding.bed` (from Stage 01)

## Planned outputs
- `outputs/panel/candidates.tsv` — heterozygome filtered to selection-eligible windows (3–10 SNPs, He ≥ 0.6)
- `outputs/panel/panel_{algo}_{size}.tsv` for `algo` in `{greedy, evenly_spaced}` and `size` in `{50, 100, 150, 200, 250}` — final panels with chrom, start, end, n_snps, He, effective_cardinality
- `reports/figures/panel_chrom_distribution.png` — Fig. 2c equivalent for Pk
- `reports/figures/panel_he_distribution.png` — per-panel He boxplot
- `reports/figures/panel_inter_marker_spacing.png` — gap-distance distribution

## Planned scripts
- `scripts/py/panel_select.py` — Python port of NB2. Argv-driven: takes heterozygome TSV + algorithm + size + chrom-allocation params, writes panel TSV.

## Exit criterion
A 100-marker panel exists with HE distribution roughly comparable to Siegel's Pv exemplar (median HE 0.70–0.81 across regions). Inter-marker spacing is reasonably uniform; no chromosome dropped entirely.

## v1 deliverable
`outputs/panel/panel_greedy_100.tsv` is the headline candidate panel. Validated in Stage 06; reported in Stage 99.
