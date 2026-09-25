# 04_window_scan.smk — sliding-window scan → the Pk heterozygome (Siegel NB1 phase 3).
#
# Slides a 200 bp window (50 bp step) along the Stage 03 discovery SNP set,
# scores every window holding ≥1 SNP, and writes the candidate-window table
# Stage 05 selects the panel from.
#
# Locked decisions L1–L4 live in `04_window_scan.README.md` and in the
# `window_scan.py` docstring. In one line each:
#   L1  emit het_raw AND het_nomiss on every row; filter on het_raw only
#   L2  het_raw on un-haploidified diploid columns; haploidify (seed 250523)
#       only for the het_nomiss branch
#   L3  frac_in_core reported, never filtered on
#   L4  windows with an identical variant-position tuple collapse, keep first
#
# Outputs are version-tagged (`window_scan.version_tag`) so a re-run lands in a
# new directory and nothing existing is overwritten. `window_scan.py` refuses to
# start if any of its three outputs already exists.
#
# Inputs:
#   outputs/qc/snps.discovery.vcf.gz (+ .csi)          (Stage 03)
#   outputs/qc/snps.summary.tsv                        (Stage 03 — count guard)
#   outputs/qc/discovery_samples.declonal.txt          (Stage 02b — cohort guard)
#   data/reference/pk_core_coding.bed                  (Stage 01 — frac_in_core)
#
# Outputs:
#   outputs/scan/<tag>/windows_all.tsv.gz              — every unique non-empty window
#   outputs/scan/<tag>/heterozygome.tsv                — ≥3 SNPs & het_raw ≥ 0.50
#   outputs/scan/<tag>/window_scan_summary.tsv         — funnel + params + per-chrom
#   reports/figures/04_window_scan/*.{pdf,png}         — 4 figures, vector + 300 dpi
#
# Runtime: minutes. 593 × 51,251 is ~60 MB as int8, so no dask (deviation D3).
# If this runs long, the cause is a per-window np.unique on an un-sliced array.

SCAN_TAG     = config["window_scan"]["version_tag"]
SCAN_OUT     = f"{SCAN_DIR}/{SCAN_TAG}"
SCAN_FIG_DIR = f"{FIG_DIR}/{config['window_scan']['figures_subdir']}"

WINDOW_SCAN_FIGURES = [
    "window_het_distribution",
    "window_het_vs_entropy",
    "window_manhattan_het",
    "window_het_raw_vs_nomiss",
]

# ---------------------------------------------------------------------------
# rule window_scan
#
# One rule, one step. Single-threaded on purpose: the work is per-chromosome
# numpy on in-memory arrays, and the machine is shared.
# ---------------------------------------------------------------------------
rule window_scan:
    input:
        vcf          = f"{QC_DIR}/snps.discovery.vcf.gz",
        vcf_index    = f"{QC_DIR}/snps.discovery.vcf.gz.csi",
        snps_summary = f"{QC_DIR}/snps.summary.tsv",
        samples      = f"{QC_DIR}/discovery_samples.declonal.txt",
        bed          = "data/reference/pk_core_coding.bed",
    output:
        windows      = f"{SCAN_OUT}/windows_all.tsv.gz",
        heterozygome = f"{SCAN_OUT}/heterozygome.tsv",
        summary      = f"{SCAN_OUT}/window_scan_summary.tsv",
    params:
        size    = config["window_scan"]["size_bp"],
        step    = config["window_scan"]["step_bp"],
        min_snps = config["window_scan"]["min_snps"],
        min_het = config["window_scan"]["min_heterozygosity"],
        seed    = config["window_scan"]["seed"],
    threads: 1
    log:
        f"{LOGS_DIR}/04_window_scan/window_scan.log",
    message:
        "04_window_scan: sliding window scan → heterozygome"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.heterozygome}) $(dirname {log})
        python scripts/py/window_scan.py \
            --vcf {input.vcf} \
            --bed {input.bed} \
            --samples {input.samples} \
            --snps-summary {input.snps_summary} \
            --window-size {params.size} \
            --step {params.step} \
            --min-snps {params.min_snps} \
            --min-het {params.min_het} \
            --seed {params.seed} \
            --out-windows {output.windows} \
            --out-heterozygome {output.heterozygome} \
            --out-summary {output.summary} \
            2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# rule plot_window_scan
#
# Four figures, each as vector PDF + 300 dpi PNG. Split from the scan so the
# figures can be re-styled without re-scanning.
# ---------------------------------------------------------------------------
rule plot_window_scan:
    input:
        windows = f"{SCAN_OUT}/windows_all.tsv.gz",
    output:
        expand(f"{SCAN_FIG_DIR}/{{fig}}.{{ext}}",
               fig=WINDOW_SCAN_FIGURES, ext=["pdf", "png"]),
    params:
        # Derived from the outputs rather than hardcoded, so the figure
        # directory always tracks the declared output paths.
        outdir   = lambda w, output: output[0].rsplit("/", 1)[0],
        min_snps = config["window_scan"]["min_snps"],
        min_het  = config["window_scan"]["min_heterozygosity"],
    threads: 1
    log:
        f"{LOGS_DIR}/04_window_scan/plot_window_scan.log",
    message:
        "04_window_scan: He distribution, He-vs-entropy, Manhattan, raw-vs-nomiss"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        python scripts/py/plot_window_scan.py \
            --windows {input.windows} \
            --min-snps {params.min_snps} \
            --min-het {params.min_het} \
            --outdir {params.outdir} \
            2>&1 | tee {log}
        """
