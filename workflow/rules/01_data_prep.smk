# 01_data_prep.smk — stage inputs and build the Pk core+coding BED.
#
# Most input staging is done by envs/install.sh (symlinks from Pop-gen, or
# user copy). This stage:
#   1. asserts every required input file is present and readable
#   2. derives the Pk core+coding BED that Stages 03 and 04 use
#   3. imports the Pop-gen Stage 2 Fws table (or recomputes if absent)
#
# Inputs: VCF (+index), reference fasta, GFF, mask list, Pop-gen Fws output.
# Outputs:
#   outputs/setup/inputs_check.ok   — sentinel file gating downstream stages
#   data/reference/pk_core_coding.bed
#   outputs/setup/fws_from_popgen.tsv

# ---------------------------------------------------------------------------
# rule check_inputs — verify every required input exists and is non-empty
# ---------------------------------------------------------------------------
rule check_inputs:
    input:
        vcf       = VCF,
        vcf_index = VCF_INDEX,
        ref_fasta = REF_FASTA,
        ref_gff   = REF_GFF,
        mask_list = MASK_LIST,
    output:
        ok = f"{SETUP_DIR}/inputs_check.ok",
    log:
        f"{LOGS_DIR}/01_data_prep/check_inputs.log",
    message:
        "01_data_prep: verifying inputs (VCF, reference, GFF, mask)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ok}) $(dirname {log})
        {{
            echo "Checking inputs..."
            for f in {input.vcf} {input.vcf_index} {input.ref_fasta} {input.ref_gff} {input.mask_list}; do
                if [[ ! -s "$f" ]]; then
                    echo "FAIL: missing or empty: $f" >&2
                    echo "Re-run envs/install.sh to symlink/copy from Pop-gen, or"
                    echo "set POPGEN_PIPELINE_PATH and re-run." >&2
                    exit 1
                fi
                echo "  OK ($(du -sh "$f" | cut -f1))  $f"
            done
            echo "All inputs present."
        }} 2>&1 | tee {log}
        touch {output.ok}
        """

# ---------------------------------------------------------------------------
# rule build_core_coding_bed
#
# Siegel's framework restricts the window scan to "core, coding regions" of
# the genome. We derive that here from the PlasmoDB GFF (CDS features) minus
# the regions in regions_to_mask.list (309 problematic regions in Pop-gen).
#
# Output is a BED file used by 03_snp_qc and 04_window_scan to restrict
# variant calls and window placement.
# ---------------------------------------------------------------------------
rule build_core_coding_bed:
    input:
        gff       = REF_GFF,
        mask_list = MASK_LIST,
        ok        = f"{SETUP_DIR}/inputs_check.ok",
    output:
        bed = "data/reference/pk_core_coding.bed",
    log:
        f"{LOGS_DIR}/01_data_prep/build_core_coding_bed.log",
    message:
        "01_data_prep: deriving Pk core+coding BED from GFF + mask list"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        bash scripts/sh/build_core_coding_bed.sh \
            {input.gff} {input.mask_list} {output.bed} \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# rule import_popgen_fws
#
# Pulls the Fws table from Pop-gen Stage 2. If POPGEN_PIPELINE_PATH is set
# and contains the file, copy it; otherwise emit a stub with instructions
# and require Stage 02 to recompute via moimix.
# ---------------------------------------------------------------------------
rule import_popgen_fws:
    input:
        ok = f"{SETUP_DIR}/inputs_check.ok",
    output:
        fws = f"{SETUP_DIR}/fws_from_popgen.tsv",
    log:
        f"{LOGS_DIR}/01_data_prep/import_popgen_fws.log",
    message:
        "01_data_prep: importing Fws table from Pop-gen Stage 2"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fws}) $(dirname {log})
        bash scripts/sh/import_popgen_fws.sh {output.fws} > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# rule audit_core_coding_bed
#
# Is `pk_core_coding.bed` actually excluding Pk's hypervariable multigene
# families? Pv4 ships a curated core-genome annotation; Pk has none, so the BED
# is built as CDS minus Pop-gen's mapping-quality mask — a mask designed for
# variant calling, not as a hypervariable-region annotation. This rule provides
# the orthogonal check: locate annotated SICAvar / kir / pseudogene loci in the
# lifted GFF and report how much of them survives into the BED and into the
# discovery SNP set.
#
# EXIT CRITERION: the script exits 2 if flagged families hold more than
# `max_flagged_snp_frac` of discovery SNPs. Non-zero fails the rule, which fails
# the workflow, which is the point — panel windows drawn from hypervariable
# families look diverse and do not amplify reliably.
#
# WHY THIS RULE IS IN STAGE 01 BUT READS A STAGE 03 OUTPUT. The audit is a
# property of the BED, which Stage 01 builds, and `DECISIONS.md` has claimed
# since 2026-08-17 that it is enforced by a Stage 01 rule. But the ">1% of
# discovery SNPs" criterion needs the discovery VCF, which Stage 03 produces.
# Snakemake resolves by DAG rather than by file order, so declaring the VCF as
# an input places the rule after Stage 03 without any circularity: nothing in
# Stages 01-03 depends on this rule's output. The rule therefore lives where
# MISSION and DECISIONS.md say it lives and still tests what matters.
#
# `audit_core_coding_bed.py --out-bed` can also emit a BED with flagged
# intervals removed. It is deliberately NOT wired: nothing consumes a cleaned
# BED, and switching Stages 03/04 onto one would change the core-genome
# definition, which is a threshold decision and out of scope here.
# ---------------------------------------------------------------------------
rule audit_core_coding_bed:
    input:
        gff = REF_GFF_LIFTED,
        bed = "data/reference/pk_core_coding.bed",
        vcf = f"{QC_DIR}/snps.discovery.vcf.gz",
    output:
        report = f"{SETUP_DIR}/core_coding_audit.tsv",
    params:
        # 1% of discovery SNPs. Hardcoded rather than read from config.yaml
        # because this run's mission does not permit editing config.yaml; the
        # value belongs under `snp_qc.max_flagged_snp_frac` and moving it there
        # is a one-line follow-up.
        max_flagged_snp_frac = 0.01,
    log:
        f"{LOGS_DIR}/01_data_prep/audit_core_coding_bed.log",
    message:
        "01_data_prep: auditing pk_core_coding.bed against SICAvar/kir/pseudogenes"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report}) $(dirname {log})
        python scripts/py/audit_core_coding_bed.py \
            --gff {input.gff} \
            --bed {input.bed} \
            --vcf {input.vcf} \
            --out-report {output.report} \
            --max-flagged-snp-frac {params.max_flagged_snp_frac} \
            2>&1 | tee {log}
        """
