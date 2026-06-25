# 03_snp_qc — SNP filtering for the discovery cohort (Siegel NB1 phase 2)

**Status:** implemented (HANDOFF #03, 2026-05-12). Not yet run.

## What this stage does
Applies Siegel et al's per-SNP filters, restricted to the monoclonal discovery cohort from Stage 02. Single bcftools pipeline (order matters: subset → filter → restrict → tag → threshold) so MAF and missingness are computed on the cohort that will actually drive panel selection.

Filters:
- biallelic SNPs only (`bcftools view -m2 -M2 -v snps`)
- `FILTER == PASS`
- restricted to `data/reference/pk_core_coding.bed` (the lifted core+coding BED from HANDOFF #02)
- MAF ≥ `snp_qc.maf_min` (default 0.10)
- F_MISSING < `snp_qc.missingness_max` (default 0.10)

## Inputs
- `data/vcf/merged_popgen.vcf.gz` (+ .csi) — Stage 01
- `outputs/qc/discovery_samples.txt` — Stage 02
- `data/reference/pk_core_coding.bed` — Stage 01

## Outputs
- `outputs/qc/snps.discovery.vcf.gz` (+ .csi) — input to Stage 04 window scan.
- `outputs/qc/snps.summary.tsv` — per-chromosome SNP counts.
- `reports/figures/snp_qc_maf_density.png`
- `reports/figures/snp_qc_missingness_density.png`

## Rules
| Rule | Purpose |
|---|---|
| `build_discovery_snp_vcf` | Single bcftools pipeline producing the final filtered VCF + CSI. |
| `summarise_snps` | Per-chromosome SNP counts + total. |
| `plot_snp_qc` | MAF + missingness densities on the kept set. |

## Exit criterion
- SNP count: ballpark 5–10K (Pk genome ~half Pv's; Siegel's Pv discovery cohort yielded ~13K).
- Per-chromosome distribution roughly flat-ish (no chromosome dominating or empty after masking).
- Density plots: MAF left-truncated at `maf_min`, F_MISSING right-truncated at `missingness_max`. Any "WARNING: N SNPs violate threshold" annotation on the plots indicates a logic bug in the filter pipeline.

## Notes
- The chromosome-naming open issue (CONTEXT-software.md #1) is resolved as of HANDOFF #02 — `pk_core_coding.bed` is now on `ordered_PKNH_*_v2` naming and `bcftools view -R` will subset cleanly.
- Intermediate VCFs are streamed (`--output-type u`) and never written to disk, so disk footprint is just the final `snps.discovery.vcf.gz` (likely <1 GiB for ~5–10K SNPs × ~615 samples).
