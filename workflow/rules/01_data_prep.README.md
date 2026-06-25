# 01_data_prep — stage inputs + build Pk core+coding BED

## What this stage does
Most input staging is done by `envs/install.sh` — it symlinks the VCF, reference, GFF, and mask from `Pop-gen_pipeline/` (or copies them if Pop-gen isn't available). This Snakemake stage does the parts that need to be re-runnable:

1. **`check_inputs`** — assert every required input exists and is non-empty. Sentinel file `outputs/setup/inputs_check.ok` gates all downstream stages.
2. **`build_core_coding_bed`** — derive `data/reference/pk_core_coding.bed` from the PlasmoDB-68 GFF (CDS features) minus the 309 problematic regions in `regions_to_mask.list`. Stages 03 and 04 use this BED to restrict variant calls and window placement to the regions Siegel's framework targets.
3. **`import_popgen_fws`** — copy the per-sample Fws table from Pop-gen Stage 2 into `outputs/setup/fws_from_popgen.tsv`. Stage 02 uses it to identify likely-monoclonal samples (Fws ≥ 0.95). If Pop-gen isn't available, the helper script emits a stub and Stage 02 falls back to recomputing via moimix.

## Inputs
| Path                                                | Source                                |
|-----------------------------------------------------|---------------------------------------|
| `data/vcf/merged_popgen.vcf.gz` (+ `.csi`)          | Pop-gen `data/vcf/`                   |
| `data/reference/PlasmoDB-67_PknowlesiA1H1_Genome.fasta` (+ indices) | Pop-gen `data/reference/PlasmoDB_version/` |
| `data/reference/PlasmoDB-68_PknowlesiA1H1.gff`      | Pop-gen `data/reference/PlasmoDB_version/` |
| `data/reference/regions_to_mask.list`               | Pop-gen `data/reference/`             |
| Pop-gen Fws table (Stage 2 output)                  | resolved by `import_popgen_fws.sh`    |

## Outputs
| Path                                          | Used by               |
|-----------------------------------------------|-----------------------|
| `outputs/setup/inputs_check.ok`               | every downstream rule |
| `data/reference/pk_core_coding.bed`           | Stage 03, Stage 04    |
| `outputs/setup/fws_from_popgen.tsv`           | Stage 02              |

## Helper scripts
- `scripts/sh/build_core_coding_bed.sh` — GFF → BED (CDS features) intersected against the inverse of the mask list. **TODO: implement.** Should use `awk` + `bedtools` for clarity.
- `scripts/sh/import_popgen_fws.sh` — locates and copies the Pop-gen Fws table; emits a clear stub error if Pop-gen isn't reachable. **TODO: implement.**

## Parameters
None at this stage — all paths come from `config.yaml`.

## Exit criterion
`outputs/setup/inputs_check.ok` exists; `data/reference/pk_core_coding.bed` exists and contains a sensible number of intervals (Pk has ~5,300 protein-coding genes; expect tens of thousands of CDS intervals after merging exons); `outputs/setup/fws_from_popgen.tsv` exists with one row per sample.

## Things that can go wrong
- **VCF index out of date** — re-run `bcftools index -c data/vcf/merged_popgen.vcf.gz` if `check_inputs` complains.
- **GFF feature names** — PlasmoDB GFF uses `CDS` for coding regions (lowercase variants exist in some annotations); helper script should be case-insensitive.
- **Mask list format** — Pop-gen's `regions_to_mask.list` uses `chrom:start-end` (one region per line). `build_core_coding_bed.sh` parses this and converts to BED.
- **Pop-gen Fws absent** — `import_popgen_fws.sh` should emit a one-line stub TSV with header only and a clear stderr message pointing the user at recomputing in Stage 02.
