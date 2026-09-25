# scripts/R/VENDORED.md

R scripts lifted from the `agnostic` pop-gen pipeline for Stage 02b. Copied,
not referenced — the pipeline must run without a sibling checkout present.

- Source tree: `Pop-gen_pipeline/agnostic`
- Commit: `e5ba5bb356d6d77388306dde0b67168ef80614dd`
- Vendored: 2026-08-28

| File | Upstream path | Role in Stage 02b |
|---|---|---|
| find_duplicates.R | scripts/R/find_duplicates.R | Replicate detection (Siegel's `Analysis_set` criterion) |
| genotype_table.R | scripts/R/genotype_table.R | VCF -> hmmIBD genotype table via contig_map |
| clonal_clusters.R | scripts/R/clonal_clusters.R | Clonal groups as connected components of the IBD graph |
| run_moimix.R | scripts/R/run_moimix.R | Fws recomputation for the A1 cross-check |
| vcf_to_gds.R | scripts/R/vcf_to_gds.R | VCF -> GDS, the input run_moimix.R needs |

All five are byte-identical to upstream apart from the vendoring header.
Re-sync by re-copying and restoring the header; do not edit in place.

## Vendored data

`data/metadata/admix_clusters.tsv` — ADMIXTURE cluster labels
(`agnostic/outputs/structure/full/admix_clusters.tsv`, same commit).
Imported as a per-sample annotation only; no pairwise IBD result is imported.
See decision L2 in the Stage 02b README.
