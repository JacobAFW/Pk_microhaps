# 02b_declonalise.smk — duplicates, IBD, and clonal thinning.
#
# Closes the last two gaps between our sample selection and Siegel's, and gives
# the discovery cohort a defensible independence assumption before any allele
# frequency is estimated from it.
#
#   1. remove replicate samples          — Siegel's third sample criterion
#                                          (`Exclusion reason == Analysis_set`),
#                                          which Stage 02 never implemented
#   2. compute IBD, derive clonal groups — hmmIBD within homogeneous clusters
#   3. optionally thin clonal groups     — one representative each, behind a flag
#
# Mirrors the `clonality.declonalize` pattern in the agnostic pop-gen pipeline:
# run everything up to IBD, detect clonality, filter to single representatives,
# re-run what requires no clones. Config vocabulary is deliberately identical
# across the two repos — same idiom in both places beats a locally neater design.
#
# Inputs:
#   data/vcf/merged_popgen.vcf.gz                 (Stage 01)
#   outputs/qc/discovery_samples.txt              (Stage 02)
#   outputs/qc/callable_fraction.tsv              (Stage 02)
#   data/metadata/admix_clusters.tsv              (vendored; decision L2)
#   data/reference/strain_preDB_version/*.fai     (contig -> integer map)
#
# Outputs:
#   outputs/qc/duplicates.remove                  — PLINK --remove list
#   outputs/qc/discovery_samples.dedup.txt        — post-dedupe cohort
#   outputs/ibd/<tag>/...                         — per-run IBD workspace
#   outputs/qc/clonal_clusters.tsv                — clonal groups
#   outputs/qc/discovery_samples.declonal.txt     — final cohort -> Stage 03
#
# See 02b_declonalise.README.md for the L1-L4 decisions and the two-run design.

import glob
import os

CLONALITY = config["clonality"]

# Run tag: outputs are keyed by the one parameter that separates Run A from
# Run B, so both can coexist on disk and be diffed directly.
#   Run A (regression) ibd_max_variant_missing: 0.10 -> outputs/ibd/geno10
#   Run B (production) ibd_max_variant_missing: 0.20 -> outputs/ibd/geno20
# Both keys accept a top-level --config override so the two runs differ by one
# flag on the command line and nothing else:
#   Run A: snakemake --config ibd_max_variant_missing=0.10
#   Run B: snakemake --config ibd_max_variant_missing=0.20
IBD_MAX_MISS = float(config.get("ibd_max_variant_missing",
                                CLONALITY["ibd_max_variant_missing"]))
DECLONALIZE = config.get("declonalize", CLONALITY["declonalize"])
if isinstance(DECLONALIZE, str):
    DECLONALIZE = DECLONALIZE.strip().lower() in ("1", "true", "yes", "on")

IBD_TAG      = "geno%02d" % round(IBD_MAX_MISS * 100)
IBD_DIR      = f"outputs/ibd/{IBD_TAG}"

# Fail loudly rather than silently producing an ADMIXTURE-shaped variant set.
# hmmIBD models IBD along the chromosome and needs linkage structure; pruning it
# away is the single easiest way to get quietly wrong answers here (L1).
LD_PRUNE = config.get("ibd_ld_prune", CLONALITY.get("ibd_ld_prune", False))
if isinstance(LD_PRUNE, str):
    LD_PRUNE = LD_PRUNE.strip().lower() in ("1", "true", "yes", "on")
if LD_PRUNE:
    raise WorkflowError(
        "clonality.ibd_ld_prune must be false — hmmIBD models IBD along the "
        "chromosome and requires linkage structure. LD pruning is an ADMIXTURE "
        "requirement, not an IBD one. See decision L1."
    )


def ibd_clusters(wildcards):
    """IBD-eligible cluster names, read from the checkpoint at DAG time."""
    ckpt = checkpoints.ibd_cluster_membership.get(**wildcards).output.keep_dir
    return sorted(
        os.path.basename(p)[: -len(".keep")]
        for p in glob.glob(os.path.join(ckpt, "*.keep"))
    )


def cluster_keep(wildcards):
    """Path to one cluster's keep-list, routed through the checkpoint.

    Taking this as a plain path string would let Snakemake build `cluster_vcf`
    without ever running the checkpoint that creates the keep-list, so asking
    for a single cluster's hmm_fract.txt on a clean tree would fail with a
    MissingInputException instead of just running the checkpoint first.
    """
    ckpt = checkpoints.ibd_cluster_membership.get().output.keep_dir
    return f"{ckpt}/{wildcards.cluster}.keep"


def ibd_fract_files(wildcards):
    return [
        f"{IBD_DIR}/{c}/{c}.hmm_fract.txt" for c in ibd_clusters(wildcards)
    ]


def ibd_fract_specs(wildcards):
    """`<path>:<cluster>` arg specs, the shape clonal_clusters.R expects."""
    return [
        f"{IBD_DIR}/{c}/{c}.hmm_fract.txt:{c}" for c in ibd_clusters(wildcards)
    ]


# ---------------------------------------------------------------------------
# rule build_contig_map
#
# hmmIBD wants an integer CHROM. genotype_table.R maps names -> integers via
# this table, replacing the hardcoded `ordered_PKNH_*_v2` gsub in the legacy
# script.
#
# NOTE: this reads the strain_A1 .fai, NOT the PlasmoDB-67 .fai that
# `reference.fasta` points at. The PlasmoDB assembly names its contigs with LT*
# accessions; the VCF uses `ordered_PKNH_*_v2`. Only the strain_A1 .fai matches
# what genotype_table.R will actually see in the VCF.
#
# Non-nuclear contigs are excluded: organellar genomes are uniparentally
# inherited and break hmmIBD's along-chromosome segment model.
# ---------------------------------------------------------------------------
rule build_contig_map:
    input:
        fai = CLONALITY["contig_fai"],
    output:
        cmap    = f"{SETUP_DIR}/contig_map.tsv",
        contigs = f"{SETUP_DIR}/nuclear_contigs.txt",
    params:
        exclude = "|".join(CLONALITY["exclude_contigs"]),
    log:
        f"{LOGS_DIR}/02b_declonalise/build_contig_map.log",
    message:
        "02b: contig -> integer map from the reference .fai (nuclear only)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cmap}) $(dirname {log})
        {{
            echo "==> Building contig map from {input.fai}"
            echo "    excluding: {params.exclude}"
            cut -f1 {input.fai} \
              | grep -Ev '^({params.exclude})$' \
              > {output.contigs}
            awk 'BEGIN{{OFS="\t"}} {{ print $1, NR }}' {output.contigs} > {output.cmap}
            echo "==> $(wc -l < {output.contigs}) nuclear contigs"
            cat {output.cmap}
        }} 2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# rule dedupe_missingness
#
# find_duplicates.R is vendored byte-identical from agnostic, so it expects
# PLINK2's .fam + .smiss. This repo has no plink2 (Stages 01-03 are pure
# bcftools), so both are synthesised from Stage 02's callable fraction:
# F_MISS = 1 - callable_fraction. See scripts/py/make_smiss.py.
# ---------------------------------------------------------------------------
rule dedupe_missingness:
    input:
        callable = f"{QC_DIR}/callable_fraction.tsv",
        samples  = f"{QC_DIR}/discovery_samples.txt",
    output:
        fam   = f"{QC_DIR}/dedupe/cohort.fam",
        smiss = f"{QC_DIR}/dedupe/cohort.smiss",
    log:
        f"{LOGS_DIR}/02b_declonalise/dedupe_missingness.log",
    message:
        "02b: per-sample missingness for the replicate tie-break"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fam}) $(dirname {log})
        python scripts/py/make_smiss.py \
            --callable {input.callable} \
            --samples {input.samples} \
            --out-fam {output.fam} \
            --out-smiss {output.smiss} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# rule remove_duplicates                                              [RULE 1]
#
# Siegel's third sample criterion. Replicates are matched on ID pattern; the
# lowest-missingness member of each group is kept and the rest written to a
# PLINK --remove list.
#
# LIMITATION (deviation A2b): matching is on ID pattern only, so this catches
# technical replicates carrying a lane suffix and nothing else. A same-patient
# re-sample under a different study code is invisible here and will surface only
# as IBD ~ 1 — indistinguishable from two patients carrying clonal parasites
# without metadata. For discovery the distinction does not matter (both should
# be thinned before frequency estimation). For the transmission interpretation
# later, it matters a great deal.
# ---------------------------------------------------------------------------
rule remove_duplicates:
    input:
        fam     = rules.dedupe_missingness.output.fam,
        smiss   = rules.dedupe_missingness.output.smiss,
        samples = f"{QC_DIR}/discovery_samples.txt",
    output:
        remove_list = f"{QC_DIR}/duplicates.remove",
        dedup  = f"{QC_DIR}/discovery_samples.dedup.txt",
    params:
        pattern = CLONALITY["duplicate_id_pattern"],
    log:
        f"{LOGS_DIR}/02b_declonalise/remove_duplicates.log",
    message:
        "02b: replicate detection (Siegel's Analysis_set criterion)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.remove_list}) $(dirname {log})
        {{
            n_in=$(wc -l < {input.samples})
            Rscript scripts/R/find_duplicates.R \
                {input.fam} {input.smiss} {output.remove_list} '{params.pattern}'

            # PLINK --remove format is "FID IID"; take IID and subtract.
            awk '{{ print $2 }}' {output.remove_list} | sort -u > {output.remove_list}.ids
            sort {input.samples} > {output.dedup}.sorted
            comm -23 {output.dedup}.sorted {output.remove_list}.ids > {output.dedup}
            rm -f {output.dedup}.sorted

            n_drop=$(wc -l < {output.remove_list})
            n_out=$(wc -l < {output.dedup})
            echo ""
            echo "==> Cohort: $n_in -> $n_out  (dropped $n_drop replicate(s))"
            echo "==> Dropped IDs:"
            cat {output.remove_list}.ids | sed 's/^/      /'

            # Stop condition: over-matching ID pattern.
            pct=$(awk -v d="$n_drop" -v t="$n_in" 'BEGIN{{ printf "%.2f", 100*d/t }}')
            echo "==> Dropped $pct% of the cohort"
            awk -v p="$pct" 'BEGIN{{ if (p > 5.0) {{
                print "ERROR: duplicate removal dropped " p "% (> 5%) — the ID "
                print "       pattern is over-matching. Stop and inspect."
                exit 1 }} }}'
        }} 2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# rule ibd_variant_set                                                [RULE 2]
#
# Decision L1: genome-wide, MAF >= 0.01, variant missingness <= configured,
# NOT LD-pruned. Three constraints drive this:
#
#   - Never the Stage 03 discovery SNPs. They are MAF >= 0.10, core-coding
#     only — a biased, marker-poor ascertainment for relatedness. Worse, they
#     are derived from the cohort we are about to thin, which makes the thinning
#     mildly circular. The IBD variant set must be independent of the
#     sample-selection step it feeds.
#   - Never LD-pruned (asserted at parse time above).
#   - More markers, but not indiscriminately: very low MAF and high-missingness
#     markers add noise rather than information.
#
# Built from the FULL VCF, subset to the post-dedupe cohort so MAF and F_MISSING
# are computed on the samples that will actually drive the analysis.
# ---------------------------------------------------------------------------
rule ibd_variant_set:
    input:
        vcf       = VCF,
        vcf_index = VCF_INDEX,
        samples   = rules.remove_duplicates.output.dedup,
        contigs   = rules.build_contig_map.output.contigs,
    output:
        vcf   = f"{IBD_DIR}/ibd_variants.vcf.gz",
        index = f"{IBD_DIR}/ibd_variants.vcf.gz.csi",
    params:
        maf      = CLONALITY["ibd_min_maf"],
        max_miss = IBD_MAX_MISS,
        min_var  = CLONALITY["ibd_min_variants"],
    threads: 4
    log:
        f"{LOGS_DIR}/02b_declonalise/ibd_variant_set_{IBD_TAG}.log",
    message:
        f"02b: IBD variant set ({IBD_TAG}) — genome-wide, no LD pruning"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        {{
            echo "==> IBD variant set [{IBD_TAG}]"
            echo "    samples:          $(wc -l < {input.samples})"
            echo "    contigs:          $(wc -l < {input.contigs}) (nuclear only)"
            echo "    MAF min:          {params.maf}"
            echo "    F_MISSING max:    {params.max_miss}"
            echo "    LD pruning:       DISABLED (required for hmmIBD)"

            regions=$(paste -sd, {input.contigs})
            bcftools view \
                  -r "$regions" \
                  --threads {threads} \
                  --output-type u \
                  {input.vcf} \
              | bcftools view \
                  --samples-file {input.samples} \
                  --force-samples \
                  --output-type u \
              | bcftools view \
                  -f PASS \
                  -m2 -M2 -v snps \
                  --output-type u \
              | bcftools +fill-tags - -- -t F_MISSING,AF \
              | bcftools view \
                  -e 'F_MISSING > {params.max_miss} || MAF < {params.maf}' \
                  --output-type z \
                  --output {output.vcf}

            bcftools index --csi --threads {threads} {output.vcf}

            n=$(bcftools index -n {output.vcf})
            echo "==> Wrote $n variants to {output.vcf}"

            # Assertion: enough markers for a meaningful IBD estimate.
            if [ "$n" -lt "{params.min_var}" ]; then
                echo "ERROR: only $n variants (< {params.min_var}). Relax "
                echo "       clonality.ibd_max_variant_missing before proceeding."
                exit 1
            fi
        }} 2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# checkpoint ibd_cluster_membership                                   [RULE 3]
#
# Decision L2. hmmIBD uses population allele frequencies, so IBD must be
# computed within genetically homogeneous groups; pooling Indonesia + Sabah into
# one population inflates apparent IBD between samples from the same region —
# population structure masquerading as relatedness.
#
# We import cluster LABELS from the agnostic pipeline, never a pairwise IBD
# result. A label is a stable per-sample annotation, and importing one keeps the
# regression comparison apples-to-apples.
#
# Samples with no label form a residual group: reported, but not clonal-called
# unless clonality.residual_clonal_calls is set. A residual group is not a
# genetically homogeneous population, so running hmmIBD on it would manufacture
# exactly the artefact L2 exists to prevent.
#
# A checkpoint because the eligible cluster list is only known after the join.
# ---------------------------------------------------------------------------
checkpoint ibd_cluster_membership:
    input:
        samples = rules.remove_duplicates.output.dedup,
        table   = CLONALITY["cluster_table"],
    output:
        assignments = f"{IBD_DIR}/cluster_assignments.tsv",
        summary     = f"{IBD_DIR}/cluster_summary.tsv",
        keep_dir    = directory(f"{IBD_DIR}/keep"),
    params:
        pattern  = CLONALITY["duplicate_id_pattern"],
        min_n    = CLONALITY["min_cluster_n"],
        residual = CLONALITY["residual_cluster_name"],
        resid_ibd = "--residual-clonal-calls" if CLONALITY.get("residual_clonal_calls") else "",
    log:
        f"{LOGS_DIR}/02b_declonalise/ibd_cluster_membership_{IBD_TAG}.log",
    message:
        "02b: cluster membership from the vendored ADMIXTURE labels (L2)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.assignments}) $(dirname {log})
        python scripts/py/assign_ibd_clusters.py \
            --samples {input.samples} \
            --cluster-table {input.table} \
            --id-pattern '{params.pattern}' \
            --min-cluster-n {params.min_n} \
            --residual-name {params.residual} \
            {params.resid_ibd} \
            --out-assignments {output.assignments} \
            --out-summary {output.summary} \
            --out-keep-dir {output.keep_dir} \
            > {log} 2>&1
        cat {log}
        """


# ---------------------------------------------------------------------------
# rule cluster_vcf
#
# Per-cluster subset with a within-cluster MAF filter. Cluster-specific MAF is
# the point of per-cluster IBD: variants polymorphic in one cluster may be fixed
# in another, and hmmIBD's allele frequencies must be the cluster's own.
# ---------------------------------------------------------------------------
rule cluster_vcf:
    input:
        vcf   = rules.ibd_variant_set.output.vcf,
        index = rules.ibd_variant_set.output.index,
        keep  = cluster_keep,
    output:
        vcf = f"{IBD_DIR}/{{cluster}}/cleaned.vcf.gz",
    params:
        maf = CLONALITY["ibd_min_maf"],
    threads: 2
    log:
        f"{LOGS_DIR}/02b_declonalise/cluster_vcf_{IBD_TAG}_{{cluster}}.log",
    message:
        "02b: per-cluster VCF + within-cluster MAF: {wildcards.cluster}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        {{
            cut -f1 {input.keep} > {output.vcf}.samples
            echo "==> Cluster {wildcards.cluster}: $(wc -l < {output.vcf}.samples) samples"
            bcftools view \
                  --samples-file {output.vcf}.samples \
                  --force-samples \
                  --threads {threads} \
                  --output-type u \
                  {input.vcf} \
              | bcftools +fill-tags - -- -t AF \
              | bcftools view \
                  -e 'MAF < {params.maf}' \
                  --output-type z \
                  --output {output.vcf}
            rm -f {output.vcf}.samples
            echo "==> $(bcftools view -H {output.vcf} | wc -l) variants after within-cluster MAF"
        }} 2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# rule genotype_table — VCF -> hmmIBD's integer-CHROM genotype table
# ---------------------------------------------------------------------------
rule genotype_table:
    input:
        vcf  = rules.cluster_vcf.output.vcf,
        cmap = rules.build_contig_map.output.cmap,
    output:
        tsv = f"{IBD_DIR}/{{cluster}}/hmmIBD_input.tsv",
    log:
        f"{LOGS_DIR}/02b_declonalise/genotype_table_{IBD_TAG}_{{cluster}}.log",
    message:
        "02b: genotype table: {wildcards.cluster}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        Rscript scripts/R/genotype_table.R {input.vcf} {input.cmap} {output.tsv} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# rule run_hmmibd                                                     [RULE 4]
#
# One thread each, mirroring agnostic. Parallelism comes from running clusters
# concurrently, not from threading a single hmmIBD call.
# ---------------------------------------------------------------------------
rule run_hmmibd:
    input:
        tsv = rules.genotype_table.output.tsv,
    output:
        fract = f"{IBD_DIR}/{{cluster}}/{{cluster}}.hmm_fract.txt",
        segs  = f"{IBD_DIR}/{{cluster}}/{{cluster}}.hmm.txt",
    params:
        prefix = f"{IBD_DIR}/{{cluster}}/{{cluster}}",
    threads: 1
    log:
        f"{LOGS_DIR}/02b_declonalise/run_hmmibd_{IBD_TAG}_{{cluster}}.log",
    message:
        "02b: hmmIBD: {wildcards.cluster}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fract}) $(dirname {log})
        hmmIBD -i {input.tsv} -o {params.prefix} > {log} 2>&1
        echo "==> $(($(wc -l < {output.fract}) - 1)) pairs" >> {log}
        """


# ---------------------------------------------------------------------------
# rule clonal_groups                                                  [RULE 5]
#
# Clonal groups are the connected components of the >= threshold pair graph,
# computed over every cluster's hmm_fract.txt at once.
# ---------------------------------------------------------------------------
rule clonal_groups:
    input:
        fracts   = ibd_fract_files,
        metadata = f"{QC_DIR}/sample_metadata.ibd.tsv",
    output:
        tsv    = f"{QC_DIR}/clonal_clusters.tsv",
        # Run-tagged copy. The canonical path above always reflects whichever
        # config last ran; this one is keyed to the parameter, so Run A and
        # Run B results coexist and can be diffed directly.
        tagged = f"{IBD_DIR}/clonal_clusters.tsv",
        focal  = (f"{IBD_DIR}/focal_{CLONALITY['focal_cluster']}_clones.tsv"
                  if CLONALITY.get("focal_cluster") else []),
    params:
        threshold = CLONALITY["clonal_ibd_threshold"],
        focal     = CLONALITY.get("focal_cluster") or "NULL",
        focal_out = (f"{IBD_DIR}/focal_{CLONALITY['focal_cluster']}_clones.tsv"
                     if CLONALITY.get("focal_cluster") else "NULL"),
        specs     = ibd_fract_specs,
    log:
        f"{LOGS_DIR}/02b_declonalise/clonal_groups_{IBD_TAG}.log",
    message:
        "02b: clonal groups as connected components of the IBD graph"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        Rscript scripts/R/clonal_clusters.R \
            {params.threshold} \
            {input.metadata} \
            {output.tsv} \
            {params.focal_out} \
            {params.focal} \
            {params.specs} \
            > {log} 2>&1
        cp {output.tsv} {output.tagged}
        cat {log}
        """


# ---------------------------------------------------------------------------
# rule ibd_metadata — clonal_clusters.R wants a `sample_id` column
# ---------------------------------------------------------------------------
rule ibd_metadata:
    input:
        meta = f"{QC_DIR}/sample_metadata.tsv",
    output:
        tsv = f"{QC_DIR}/sample_metadata.ibd.tsv",
    log:
        f"{LOGS_DIR}/02b_declonalise/ibd_metadata.log",
    message:
        "02b: renaming sample -> sample_id for clonal_clusters.R"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        awk 'BEGIN{{FS=OFS="\t"}}
             NR==1 {{ for(i=1;i<=NF;i++){{ if($i=="sample") $i="sample_id";
                                          if($i=="region_geo") $i="geography" }} }}
             {{ print }}' {input.meta} > {output.tsv}
        echo "wrote {output.tsv} ($(($(wc -l < {output.tsv}) - 1)) samples)" | tee {log}
        """


# ---------------------------------------------------------------------------
# rule declonalised_cohort                                            [RULE 6]
#
# Decision L3. One representative per clonal group, chosen by lowest missingness
# — the same tie-break find_duplicates.R uses for technical replicates.
#
# `declonalize: false` by default so the headline v1 numbers stay comparable to
# what is already documented; the flag makes the sensitivity analysis one
# parameter away. When false this is a pass-through copy, so the DAG shape is
# identical either way.
# ---------------------------------------------------------------------------
rule declonalised_cohort:
    input:
        samples = rules.remove_duplicates.output.dedup,
        clonal  = rules.clonal_groups.output.tsv,
        smiss   = rules.dedupe_missingness.output.smiss,
    output:
        samples = f"{QC_DIR}/discovery_samples.declonal.txt",
        report  = f"{QC_DIR}/declonalisation_report.tsv",
    params:
        flag = "--declonalize" if DECLONALIZE else "",
    log:
        f"{LOGS_DIR}/02b_declonalise/declonalised_cohort.log",
    message:
        "02b: applying the thinning policy (L3)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.samples}) $(dirname {log})
        python scripts/py/apply_declonalisation.py \
            --samples {input.samples} \
            --clonal-clusters {input.clonal} \
            --smiss {input.smiss} \
            {params.flag} \
            --out-samples {output.samples} \
            --out-report {output.report} \
            > {log} 2>&1
        cat {log}
        """


# ---------------------------------------------------------------------------
# rule plot_ibd_distribution
#
# The fract_sites_IBD histogram is the health check on the whole stage: bulk
# near 0, small spike near 1. Substantial mass in the middle means residual
# polyclonality or a population-structure artefact from L2 — stop and inspect.
# ---------------------------------------------------------------------------
rule plot_ibd_distribution:
    input:
        fracts = ibd_fract_files,
    output:
        png = f"{FIG_DIR}/ibd_fract_distribution_{IBD_TAG}.png",
    params:
        threshold = CLONALITY["clonal_ibd_threshold"],
        tag       = IBD_TAG,
    log:
        f"{LOGS_DIR}/02b_declonalise/plot_ibd_distribution_{IBD_TAG}.log",
    message:
        "02b: fract_sites_IBD distribution"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.png}) $(dirname {log})
        python scripts/py/plot_ibd.py \
            --fract {input.fracts} \
            --threshold {params.threshold} \
            --tag {params.tag} \
            --out {output.png} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# rule fws_cross_check                                                [RULE 7]
#
# Deviation A1. Our Fws is imported from Pop-gen's `Proportion` column and has
# never been verified. run_moimix.R computes Fws properly; compare_fws_sources.py
# reports Pearson r, max absolute difference, and — the actual verdict — how many
# samples cross the 0.95 threshold under one source but not the other. If that
# number is zero, `Proportion` is fit for purpose whatever it is called.
#
# NOT in FINAL_TARGETS. Status as of 2026-08-28: SeqArray, gdsfmt and BiocParallel
# were installed successfully, so `vcf_to_gds` runs. `moimix` is still blocked —
# it hard-depends on SeqVarTools, whose chain
# (SeqVarTools -> logistf -> mice -> lme4 -> RcppEigen/nloptr) will not build
# against R 4.3.1 here: RcppEigen fails with "no member named 'Rlog1p' in
# namespace 'std'", an R-headers macro collision that -DR_NO_REMAP did not fix.
# Bioconductor-source, conda (no osx-arm64 build), GitHub-upstream and
# dependencies=FALSE routes were all tried.
#
# Build it explicitly once moimix is available:
#     snakemake outputs/qc/fws_cross_check.tsv
# ---------------------------------------------------------------------------
rule vcf_to_gds:
    input:
        vcf     = rules.ibd_variant_set.output.vcf,
        samples = rules.remove_duplicates.output.dedup,
    output:
        gds = f"{IBD_DIR}/fws/cohort.gds",
    log:
        f"{LOGS_DIR}/02b_declonalise/vcf_to_gds.log",
    message:
        "02b: VCF -> GDS for moimix"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gds}) $(dirname {log})
        Rscript scripts/R/vcf_to_gds.R {input.vcf} {output.gds} > {log} 2>&1
        """


rule fws_cross_check:
    input:
        gds      = rules.vcf_to_gds.output.gds,
        imported = f"{SETUP_DIR}/fws_from_popgen.tsv",
    output:
        fws    = f"{QC_DIR}/fws_computed.tsv",
        baf    = f"{QC_DIR}/fws_baf.tsv",
        report = f"{QC_DIR}/fws_cross_check.tsv",
    params:
        threshold = config["sample_qc"]["fws_threshold"],
        seed      = 42,
    log:
        f"{LOGS_DIR}/02b_declonalise/fws_cross_check.log",
    message:
        "02b: Fws cross-check — the A1 verdict"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report}) $(dirname {log})
        {{
            Rscript scripts/R/run_moimix.R {input.gds} {output.fws} {output.baf} {params.seed}
            python scripts/py/compare_fws_sources.py \
                --imported {input.imported} \
                --computed {output.fws} \
                --threshold {params.threshold} \
                --out-report {output.report} \
                --fail-on-crossing
        }} 2>&1 | tee {log}
        """
