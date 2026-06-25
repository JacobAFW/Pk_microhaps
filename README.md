# P. knowlesi microhaplotype panel discovery

A reproducible Snakemake pipeline for discovering a candidate microhaplotype
panel for *Plasmodium knowlesi* (Pk) from whole-genome sequencing data.

## What this repo is

The pipeline ports the *P. vivax* microhaplotype discovery framework of
[Siegel et al. 2024 (*Nat Commun*)](https://www.nature.com/articles/s41467-024-51015-3)
onto Pk WGS. Each phase of the source framework maps to a Snakemake stage;
stage outputs feed forward and a Quarto report grows one chapter per completed
stage. The year-one deliverable is a candidate ~100-microhaplotype panel plus a
[`paneljudge`](https://github.com/aimeertaylor/paneljudge) IBD-RMSE validation.
Stack: Snakemake, `bcftools`/`samtools`, Python (`scikit-allel`), and R
(`paneljudge`, `moimix`). Stages 01–03 are implemented; 04–06 are scaffolded.

## What's included — and what's deliberately not

This repository contains **code only**. By design it does **not** include:

- raw or processed **data** (sequence, genotype, phenotype, tabular records) —
  VCFs, reference genomes/indices, BEDs, GFFs, and all pipeline outputs;
- **sample sheets, manifests, or metadata** that link samples to individuals;
- any **identifying or sensitive** information — in particular the entire
  sample-prioritisation workstream, which operates on patient mastersheets and
  case-level epidemiological data (including geographic coordinates) and whose
  scripts embed patient-level filter logic, named individuals, and budget
  figures;
- **internal project notes** — drafts, handoffs, and status/planning/context/
  memory documents containing collaborator names and unpublished detail;
- **credentials, tokens, or environment files**, and machine-generated local-path
  files (e.g. `envs/activate.sh`);
- vendored toolchains/caches and third-party cloned repos / copyrighted papers.

Data lives outside version control (institutional storage / controlled-access
source); the scripts expect it at the paths described in `workflow/config.yaml`.

## Reproducing the analysis

1. **Environment:** `bash envs/install.sh` bootstraps a project-scoped
   conda/mamba environment (nothing is written outside the project directory).
   A pinned `envs/environment.lock.yaml` (osx-arm64) records exact versions.
2. **Inputs:** provide your own under `data/` (paths configurable in
   `workflow/config.yaml`): a multi-sample Pk VCF (`+ .csi`), a PlasmoDB Pk
   reference FASTA (`+ indices`) and GFF, a region-mask list, and a per-sample
   Fws table. None are distributed here.
3. **Run:** `snakemake -n` to preview the DAG, then `snakemake --cores N`.
   All stage thresholds live in `workflow/config.yaml`.

## Structure

```
workflow/
  Snakefile              # top-level workflow
  config.yaml            # parameters and (relative) paths
  rules/*.smk            # per-stage Snakemake rules (01, 02, 03, 99)
  rules/*.README.md      # per-stage documentation (01–06, 99)
scripts/
  py/                    # Python stage scripts (scikit-allel based)
  sh/                    # shell helpers (BED construction, imports)
envs/
  install.sh             # one-shot environment bootstrap
  environment.lock.yaml  # pinned conda environment (osx-arm64)
_quarto.yml              # Quarto report configuration
```

## License

<choose one — e.g. MIT — before publishing>

## Citation

The discovery and selection methodology follows Siegel et al. 2024
(*Nat Commun*); panel evaluation uses `paneljudge`. Please cite the original
work when using this pipeline.
