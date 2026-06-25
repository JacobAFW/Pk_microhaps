# 04_window_scan — Pk heterozygome via sliding window scan (Siegel NB1 phase 3)

**Status:** stub — rules not yet written.

## What this stage will do
Reproduce Siegel's `1_evaluating_marker_candidates.ipynb` for Pk. Slide 200 bp windows (50 bp step) across the core+coding BED, score each candidate window, filter to the heterozygome.

## Planned inputs
- `outputs/qc/snps.discovery.vcf.gz` (from Stage 03)
- `data/reference/pk_core_coding.bed` (from Stage 01)
- `outputs/qc/discovery_samples.txt` (from Stage 02)

## Planned outputs
- `outputs/scan/all_windows.tsv` — every 200 bp window with ≥1 SNP, columns: chrom, start, end, n_snps, mean_maf, heterozygosity, effective_cardinality
- `outputs/scan/heterozygome.tsv` — the filtered subset (≥3 SNPs AND He ≥ 0.5)
- `reports/figures/heterozygome_manhattan.png` — Fig. 1c equivalent for Pk
- `reports/figures/heterozygome_per_chrom_count.png` — sanity bar plot

## Planned scripts
- `scripts/py/window_scan.py` — Python port of NB1 phase 3 using `scikit-allel`. Takes the VCF + sample list + BED + window params from `config.yaml`. Argv-driven: callable from Snakemake or standalone.

## Exit criterion
Heterozygome count is in the same order of magnitude as Siegel's Pv numbers (3,830 windows on Pv4; for Pk we'd expect roughly proportional to genome size and SNP density — say 1,000–5,000 windows). Manhattan plot shows windows distributed across all 14 chromosomes with the expected end-of-chromosome density bias Siegel noted.

## Reference implementation
- `external/clones/vivax-mhaps/notebooks/1_evaluating_marker_candidates.ipynb` — the original notebook to port. Uses `scikit-allel` for VCF traversal and per-window stats.
