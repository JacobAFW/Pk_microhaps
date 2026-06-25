# 02_sample_qc.smk — monoclonal discovery-cohort selection (Siegel NB1 phase 1).
#
# Filter the 990-sample Pk WGS cohort to high-quality monoclonal samples for
# microhap discovery. Two criteria, both Siegel's:
#   Fws ≥ 0.95                            (monoclonal proxy; from Pop-gen Stage 2)
#   callable_fraction ≥ 0.50              (≥50% of pk_core_coding.bed SNPs called)
#
# Inputs:
#   data/vcf/merged_popgen.vcf.gz                    (Stage 01)
#   data/reference/pk_core_coding.bed                (Stage 01)
#   outputs/setup/fws_from_popgen.tsv                (Stage 01)
#   data/metadata/samples.tsv                        (this stage — imported from Pop-gen)
#
# Outputs:
#   outputs/qc/callable_fraction.tsv                 (sample, n_called, n_total, callable_fraction)
#   outputs/qc/sample_metadata.tsv                   (sample, fws, callable_fraction, region_geo, include_for_discovery)
#   outputs/qc/discovery_samples.txt                 (one sample ID per line)
#   reports/figures/fws_distribution.png
#   reports/figures/callable_fraction_distribution.png
#
# Open caveats (logged in CONTEXT-software.md → Open issues #2):
#   - "Proportion" column in fws_from_popgen.tsv is assumed equivalent to Fws.
#     merge_sample_qc.py asserts 0 ≤ values ≤ 1; flag a TODO for confirmation.
#   - callable_fraction here is a proxy: fraction of SNP sites in core+coding
#     with non-missing genotype, not coverage-based. Conservative when paired
#     with Fws ≥ 0.95. TODO: revisit with BAM-derived callable if a coverage
#     pipeline lands.

# ---------------------------------------------------------------------------
# rule import_popgen_metadata — pull sample metadata (country/state) from Pop-gen
# ---------------------------------------------------------------------------
rule import_popgen_metadata:
    input:
        ok = f"{SETUP_DIR}/inputs_check.ok",
    output:
        meta = "data/metadata/samples.tsv",
    log:
        f"{LOGS_DIR}/02_sample_qc/import_popgen_metadata.log",
    message:
        "02_sample_qc: importing sample metadata (country/state) from Pop-gen"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.meta}) $(dirname {log})
        bash scripts/sh/import_popgen_metadata.sh {output.meta} > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# rule compute_callable_fraction
#
# Per-sample fraction of SNPs in pk_core_coding.bed with non-missing genotype.
# This is the v1 proxy for Siegel's coverage-based callable; tight enough when
# combined with Fws ≥ 0.95 (low-coverage samples don't pass moimix's Fws calc).
#
# Uses bcftools query for streaming + awk for the per-sample tally. No
# intermediate VCFs.
# ---------------------------------------------------------------------------
rule compute_callable_fraction:
    input:
        vcf       = VCF,
        vcf_index = VCF_INDEX,
        bed       = "data/reference/pk_core_coding.bed",
    output:
        tsv = f"{QC_DIR}/callable_fraction.tsv",
    log:
        f"{LOGS_DIR}/02_sample_qc/compute_callable_fraction.log",
    threads: 2
    message:
        "02_sample_qc: per-sample callable fraction over pk_core_coding.bed"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        {{
            echo "==> Computing per-sample callable fraction over $(wc -l < {input.bed}) regions"
            # bcftools query streams [SAMPLE\tGT] per genotype; awk tallies per
            # sample. Missing genotypes are ./. or . — count anything else as called.
            bcftools query -R {input.bed} -f '[%SAMPLE\t%GT\n]' {input.vcf} \
              | awk 'BEGIN{{OFS="\t"}}
                     {{
                         tot[$1]++
                         if ($2 != "./." && $2 != "." && $2 != ".|.") called[$1]++
                     }}
                     END {{
                         print "sample","n_called","n_total","callable_fraction"
                         for (s in tot) {{
                             cf = (tot[s] > 0) ? called[s] / tot[s] : 0
                             printf "%s\t%d\t%d\t%.6f\n", s, called[s]+0, tot[s], cf
                         }}
                     }}' \
              > {output.tsv}
            n=$(($(wc -l < {output.tsv}) - 1))
            echo "==> Wrote $n samples to {output.tsv}"
        }} 2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# rule merge_sample_qc
#
# Join Fws + callable + metadata, apply both thresholds, write the discovery
# sample list. All sanity assertions live in the python script (Fws range,
# sample-ID overlap with VCF, region_geo populated).
# ---------------------------------------------------------------------------
rule merge_sample_qc:
    input:
        fws      = f"{SETUP_DIR}/fws_from_popgen.tsv",
        callable = f"{QC_DIR}/callable_fraction.tsv",
        meta     = "data/metadata/samples.tsv",
    output:
        sample_meta = f"{QC_DIR}/sample_metadata.tsv",
        discovery   = f"{QC_DIR}/discovery_samples.txt",
    params:
        fws_threshold      = config["sample_qc"]["fws_threshold"],
        callable_threshold = config["sample_qc"]["callable_threshold"],
    log:
        f"{LOGS_DIR}/02_sample_qc/merge_sample_qc.log",
    message:
        "02_sample_qc: merging Fws + callable + metadata; applying thresholds"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.sample_meta}) $(dirname {log})
        python scripts/py/merge_sample_qc.py \
            --fws {input.fws} \
            --callable {input.callable} \
            --metadata {input.meta} \
            --fws-threshold {params.fws_threshold} \
            --callable-threshold {params.callable_threshold} \
            --out-metadata {output.sample_meta} \
            --out-discovery {output.discovery} \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# rule plot_sample_qc — Fws + callable distributions for the sanity-check report
# ---------------------------------------------------------------------------
rule plot_sample_qc:
    input:
        sample_meta = f"{QC_DIR}/sample_metadata.tsv",
    output:
        fws_plot      = f"{FIG_DIR}/fws_distribution.png",
        callable_plot = f"{FIG_DIR}/callable_fraction_distribution.png",
    params:
        fws_threshold      = config["sample_qc"]["fws_threshold"],
        callable_threshold = config["sample_qc"]["callable_threshold"],
    log:
        f"{LOGS_DIR}/02_sample_qc/plot_sample_qc.log",
    message:
        "02_sample_qc: plotting Fws + callable distributions"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fws_plot}) $(dirname {log})
        python scripts/py/plot_sample_qc.py \
            --sample-metadata {input.sample_meta} \
            --fws-threshold {params.fws_threshold} \
            --callable-threshold {params.callable_threshold} \
            --out-fws {output.fws_plot} \
            --out-callable {output.callable_plot} \
            > {log} 2>&1
        """
