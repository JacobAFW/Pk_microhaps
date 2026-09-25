# 03_snp_qc.smk — SNP-level QC for the discovery cohort (Siegel NB1 phase 2).
#
# Apply Siegel et al's per-SNP filters, restricted to the monoclonal sample
# set from Stage 02:
#   FILTER == PASS                         (caller's hard filter)
#   biallelic SNPs                          (Pk panel is biallelic-SNP only)
#   MAF ≥ 0.10  (config: snp_qc.maf_min)
#   F_MISSING < 0.10  (config: snp_qc.missingness_max)
#   intervals ∩ pk_core_coding.bed
#
# The bcftools pipeline is run as a single rule (subset → filter → restrict →
# MAF/missingness) because the intermediate VCFs aren't independently useful
# scientific artefacts. Summarisation + plotting are split out as separate
# rules so they can be re-run cheaply.
#
# Inputs:
#   data/vcf/merged_popgen.vcf.gz                     (Stage 01)
#   outputs/qc/discovery_samples.declonal.txt         (Stage 02b)
#   data/reference/pk_core_coding.bed                 (Stage 01)
#
# Outputs:
#   outputs/qc/snps.discovery.vcf.gz  (+ .csi)        — input to Stage 04
#   outputs/qc/snps.summary.tsv                       — per-chrom SNP counts
#   reports/figures/snp_qc_maf_density.png
#   reports/figures/snp_qc_missingness_density.png

# ---------------------------------------------------------------------------
# rule build_discovery_snp_vcf
#
# Single bcftools pipeline. Order matters:
#   1. subset to discovery samples (so MAF is computed on the cohort that
#      will actually drive panel selection). Since Stage 02b, that list is the
#      post-dedupe (and optionally declonalised) cohort — MAF and F_MISSING
#      recompute automatically from the dependency edge, so flipping
#      clonality.declonalize re-derives the SNP set with no orchestration.
#   2. PASS + biallelic SNPs
#   3. restrict to pk_core_coding.bed
#   4. fill F_MISSING and AF tags on the current sample set
#   5. apply MAF + missingness thresholds
# ---------------------------------------------------------------------------
rule build_discovery_snp_vcf:
    input:
        vcf       = VCF,
        vcf_index = VCF_INDEX,
        bed       = "data/reference/pk_core_coding.bed",
        samples   = f"{QC_DIR}/discovery_samples.declonal.txt",
    output:
        vcf   = f"{QC_DIR}/snps.discovery.vcf.gz",
        index = f"{QC_DIR}/snps.discovery.vcf.gz.csi",
    params:
        maf_min         = config["snp_qc"]["maf_min"],
        missingness_max = config["snp_qc"]["missingness_max"],
    threads: 4
    log:
        f"{LOGS_DIR}/03_snp_qc/build_discovery_snp_vcf.log",
    message:
        "03_snp_qc: discovery VCF — PASS + biallelic + MAF + missingness + core+coding"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        {{
            echo "==> Building discovery SNP VCF"
            echo "    samples:     $(wc -l < {input.samples})"
            echo "    intervals:   $(wc -l < {input.bed})"
            echo "    MAF min:     {params.maf_min}"
            echo "    F_MISSING <  {params.missingness_max}"

            # NOTE: -R requires a random-access (indexed) source, so the BED
            # restriction must run against the original VCF — it can't read
            # the BED restriction from stdin. Sample subset comes second so
            # MAF/F_MISSING are still computed on the discovery cohort.
            bcftools view \
                  -R {input.bed} \
                  --threads {threads} \
                  --output-type u \
                  {input.vcf} \
              | bcftools view \
                  --samples-file {input.samples} \
                  --output-type u \
              | bcftools view \
                  -f PASS \
                  -m2 -M2 -v snps \
                  --output-type u \
              | bcftools +fill-tags - -- -t F_MISSING,AF \
              | bcftools view \
                  -e 'F_MISSING > {params.missingness_max} || MAF < {params.maf_min}' \
                  --output-type z \
                  --output {output.vcf}

            bcftools index --csi --threads {threads} {output.vcf}

            n=$(bcftools view -H {output.vcf} | wc -l)
            echo "==> Wrote $n SNPs to {output.vcf}"
        }} 2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# rule summarise_snps
#
# Per-chromosome SNP count + total. Cheap sanity check before window scanning.
# ---------------------------------------------------------------------------
rule summarise_snps:
    input:
        vcf   = f"{QC_DIR}/snps.discovery.vcf.gz",
        index = f"{QC_DIR}/snps.discovery.vcf.gz.csi",
    output:
        tsv = f"{QC_DIR}/snps.summary.tsv",
    log:
        f"{LOGS_DIR}/03_snp_qc/summarise_snps.log",
    message:
        "03_snp_qc: per-chromosome SNP counts"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        {{
            echo "==> Summarising SNP distribution by chromosome"
            {{
                printf "chrom\tn_snps\n"
                bcftools query -f '%CHROM\n' {input.vcf} \
                  | sort \
                  | uniq -c \
                  | awk 'BEGIN{{OFS="\t"}} {{ print $2, $1 }}'
                total=$(bcftools view -H {input.vcf} | wc -l)
                printf "TOTAL\t%d\n" "$total"
            }} > {output.tsv}
            echo "==> Wrote {output.tsv}"
            cat {output.tsv}
        }} 2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# rule plot_snp_qc
#
# MAF + per-SNP missingness densities on the kept set. Confirms the filters
# bit and the kept distribution looks sane (MAF should be left-truncated at
# 0.10, missingness right-truncated at 0.10).
# ---------------------------------------------------------------------------
rule plot_snp_qc:
    input:
        vcf   = f"{QC_DIR}/snps.discovery.vcf.gz",
        index = f"{QC_DIR}/snps.discovery.vcf.gz.csi",
    output:
        maf_plot     = f"{FIG_DIR}/snp_qc_maf_density.png",
        miss_plot    = f"{FIG_DIR}/snp_qc_missingness_density.png",
    params:
        maf_min         = config["snp_qc"]["maf_min"],
        missingness_max = config["snp_qc"]["missingness_max"],
    log:
        f"{LOGS_DIR}/03_snp_qc/plot_snp_qc.log",
    message:
        "03_snp_qc: MAF + missingness density plots on kept SNPs"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.maf_plot}) $(dirname {log})
        python scripts/py/plot_snp_qc.py \
            --vcf {input.vcf} \
            --maf-min {params.maf_min} \
            --missingness-max {params.missingness_max} \
            --out-maf {output.maf_plot} \
            --out-missingness {output.miss_plot} \
            > {log} 2>&1
        """
