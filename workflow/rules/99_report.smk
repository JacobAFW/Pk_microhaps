# 99_report.smk — render the project's Quarto book.
#
# Book config: _quarto.yml (project root).
# Chapters:    reports/qmd/*.qmd  (one per completed pipeline stage).
# Output:      reports/_book/index.html (multi-page HTML site, sidebar navigation).
#
# The rule depends on every output that the chapters embed. As stages are
# added, extend the `input:` block with the new stage's outputs + chapter
# source so Snakemake correctly re-renders when any of them change.
#
# Currently in scope: Stages 01, 02, 03. Chapters for not-yet-implemented
# stages (04+) are deliberately absent from _quarto.yml; add them as the
# stages land.

rule render_report:
    input:
        # ---- Quarto sources -----------------------------------------------
        config   = "_quarto.yml",
        index    = "index.qmd",
        ch01     = "reports/qmd/01_data_prep.qmd",
        ch02     = "reports/qmd/02_sample_qc.qmd",
        ch03     = "reports/qmd/03_snp_qc.qmd",
        # ---- Stage 01 outputs the chapters read --------------------------
        bed              = "data/reference/pk_core_coding.bed",
        fws              = f"{SETUP_DIR}/fws_from_popgen.tsv",
        inputs_ok        = f"{SETUP_DIR}/inputs_check.ok",
        # ---- Stage 02 outputs the chapters read --------------------------
        sample_meta      = f"{QC_DIR}/sample_metadata.tsv",
        discovery        = f"{QC_DIR}/discovery_samples.txt",
        callable_tsv     = f"{QC_DIR}/callable_fraction.tsv",
        fws_plot         = f"{FIG_DIR}/fws_distribution.png",
        callable_plot    = f"{FIG_DIR}/callable_fraction_distribution.png",
        # ---- Stage 03 outputs the chapters read --------------------------
        snps_vcf         = f"{QC_DIR}/snps.discovery.vcf.gz",
        snps_summary     = f"{QC_DIR}/snps.summary.tsv",
        maf_plot         = f"{FIG_DIR}/snp_qc_maf_density.png",
        miss_plot        = f"{FIG_DIR}/snp_qc_missingness_density.png",
    output:
        # Quarto book entry point. Treated as the canonical render-fresh marker.
        html = "reports/_book/index.html",
    log:
        f"{LOGS_DIR}/99_report/render_report.log",
    message:
        "99_report: rendering Quarto book → reports/_book/"
    shell:
        r"""
        set -euo pipefail
        mkdir -p reports/_book $(dirname {log})
        {{
            echo "==> Rendering Quarto book"
            quarto --version
            # Build from project root so execute-dir: project resolves correctly.
            quarto render
            echo "==> Wrote reports/_book/index.html"
        }} 2>&1 | tee {log}
        """
