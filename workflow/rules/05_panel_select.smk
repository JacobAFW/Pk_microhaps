# 05_panel_select.smk — candidate panel selection (Siegel NB2).
#
# **STATUS: authored and linted 2026-09-25, NEVER RUN.** The include is live and
# the DAG resolves, but no panel has been produced or reviewed. Exit criteria in
# `05_panel_select.README.md` are untested.
#
# Takes the Stage 04 heterozygome, applies the tighter selection filters
# (3-10 SNPs, He >= 0.60), collapses overlapping windows into independent loci,
# allocates markers across chromosomes by length, and emits panels of each
# configured size under both algorithms.
#
# The load-bearing detail: Stage 04 windows are 200 bp at a 50 bp step, so
# neighbouring candidates overlap and describe the same genome. On the Stage 04
# output, 19,674 selection-eligible windows collapse to 5,915 independent loci.
# Selection runs over loci and enforces `panel_select.min_spacing_bp`, so a panel
# cannot spend four markers on one locus.
#
# Inputs:
#   outputs/scan/<scan_tag>/heterozygome.tsv                 (Stage 04)
#   data/reference/strain_preDB_version/...Icor.fasta.fai    (contig lengths)
#
# Outputs (version-tagged; nothing existing is overwritten):
#   outputs/panel/<tag>/candidates.tsv
#   outputs/panel/<tag>/panel_{algorithm}_{size}.tsv
#
# The v1 deliverable is `panel_greedy_100.tsv`.

PANEL_TAG  = config["panel_select"]["version_tag"]
PANEL_OUT  = f"{PANEL_DIR}/{PANEL_TAG}"
PANEL_ALGOS = config["panel_select"]["algorithms"]
PANEL_SIZES = config["panel_select"]["panel_sizes"]

# ---------------------------------------------------------------------------
# rule select_panel
#
# One rule per (algorithm, size). Each is seconds of work on a ~6k-row table, so
# there is no reason to batch them; separate outputs keep the DAG legible and let
# a single panel be rebuilt without touching the others.
# ---------------------------------------------------------------------------
rule select_panel:
    input:
        heterozygome = f"{SCAN_OUT}/heterozygome.tsv",
        fai          = config["clonality"]["contig_fai"],
    output:
        candidates = f"{PANEL_OUT}/candidates.{{algorithm}}_{{size}}.tsv",
        panel      = f"{PANEL_OUT}/panel_{{algorithm}}_{{size}}.tsv",
    params:
        min_snps     = config["panel_select"]["snp_min_per_window"],
        max_snps     = config["panel_select"]["snp_max_per_window"],
        min_het      = config["panel_select"]["heterozygosity_min"],
        min_spacing  = config["panel_select"]["min_spacing_bp"],
    wildcard_constraints:
        algorithm = "|".join(PANEL_ALGOS),
        size      = r"\d+",
    threads: 1
    log:
        f"{LOGS_DIR}/05_panel_select/select_panel_{{algorithm}}_{{size}}.log",
    message:
        "05_panel_select: {wildcards.algorithm} panel of {wildcards.size} markers"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.panel}) $(dirname {log})
        python scripts/py/panel_select.py \
            --heterozygome {input.heterozygome} \
            --fai {input.fai} \
            --min-snps {params.min_snps} \
            --max-snps {params.max_snps} \
            --min-het {params.min_het} \
            --min-spacing-bp {params.min_spacing} \
            --algorithm {wildcards.algorithm} \
            --panel-size {wildcards.size} \
            --out-candidates {output.candidates} \
            --out-panel {output.panel} \
            2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# rule panel_select_all
#
# Convenience aggregate: every algorithm x size combination. Not in
# FINAL_TARGETS — Stage 05 is scaffolded, not commissioned, so `snakemake` with
# no target must not start building panels. Run it explicitly:
#     snakemake --cores 3 panel_select_all
# ---------------------------------------------------------------------------
rule panel_select_all:
    input:
        expand(f"{PANEL_OUT}/panel_{{algorithm}}_{{size}}.tsv",
               algorithm=PANEL_ALGOS, size=PANEL_SIZES),
    message:
        "05_panel_select: all {} algorithm x size combinations".format(
            len(PANEL_ALGOS) * len(PANEL_SIZES))
