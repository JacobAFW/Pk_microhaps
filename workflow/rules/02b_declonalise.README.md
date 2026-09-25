# Stage 02b — duplicates, IBD, and clonal thinning

Closes the last two gaps between our sample selection and Siegel's, and gives the
discovery cohort a defensible independence assumption before any allele frequency
is estimated from it.

Three things, in order:

1. **Remove replicate samples** — Siegel's third sample criterion, which Stage 02
   never implemented.
2. **Compute IBD and derive clonal groups** — hmmIBD within genetically
   homogeneous clusters.
3. **Optionally thin each clonal group** to one representative, behind a config
   flag, letting Stage 03+ re-run on the thinned cohort.

This mirrors the `clonality.declonalize` pattern in the `agnostic` pop-gen
pipeline: run everything up to IBD, detect clonality, filter to single
representatives, re-run what requires no clones. **The config vocabulary is
deliberately identical in both repos** — the same idiom in both places is worth
more than a locally neater design.

---

## Why this is worth a stage of its own

- **Siegel's sample criteria were three, not two.** Fws, callable, *and*
  `Exclusion reason == Analysis_set`, which drops technical and biological
  duplicates. Stage 02 implemented two. See `CONTEXT-method-deviations.md` A2.
- **Clonal over-representation distorts window ranking.** He = 1 − Σf², so a
  repeated haplotype pushes Σf² up and He down. The direction is conservative —
  you miss good windows rather than selecting bad ones — but the deflation
  concentrates in windows where the repeated genotype is distinctive, so it
  distorts the *ranking*. Panel selection is a ranking exercise (top-N per
  region), so this produces a suboptimal panel, not a wrong one. It also acts
  upstream at Stage 03's MAF ≥ 0.10 gate, changing which SNPs enter at all.
- **It was meant to close deviation A1 as a side-effect** — see rule 7. The
  cross-check is written and wired but **has not run**: `moimix` will not install
  in this environment (`SeqVarTools → logistf → mice → lme4 → RcppEigen`, which
  fails under R 4.3.1 with `no member named 'Rlog1p' in namespace 'std'`).
  A1 remains open, narrowed to that one install.

---

## Provenance

Four R scripts and one metadata table are **vendored** from the `agnostic`
pop-gen pipeline at commit `e5ba5bb356d6d77388306dde0b67168ef80614dd`. They are
byte-identical to upstream apart from a provenance header. See
`scripts/R/VENDORED.md`.

They are *copied*, not referenced: the pipeline must run without a sibling
checkout, which a collaborator will not have.

| Vendored | Role |
|---|---|
| `scripts/R/find_duplicates.R` | Replicate detection — Siegel's `Analysis_set` filter |
| `scripts/R/genotype_table.R` | VCF → hmmIBD table via a contig_map |
| `scripts/R/clonal_clusters.R` | Clonal groups as connected components of the IBD graph |
| `scripts/R/run_moimix.R` | Fws recomputation (A1 cross-check) |
| `scripts/R/vcf_to_gds.R` | VCF → GDS, the input `run_moimix.R` needs |
| `data/metadata/admix_clusters.tsv` | ADMIXTURE cluster labels (decision L2) |

`agnostic`'s `ibd.clonal_ibd_threshold` is already `0.95` — exactly Siegel's
clone threshold. Nothing needed reconciling; the number is kept and both are
cited.

---

## Decisions

### L1 — Which variants for IBD
**Genome-wide, MAF ≥ 0.01, variant missingness ≤ `ibd_max_variant_missing`, NOT
LD-pruned.** Three constraints:

- **Never the Stage 03 discovery SNPs.** They are MAF ≥ 0.10, core-coding only —
  a biased, marker-poor ascertainment for relatedness. Worse, they are derived
  *from the cohort we are about to thin*, which makes the thinning mildly
  circular. The IBD variant set must be independent of the sample-selection step
  it feeds.
- **Never LD-pruned.** hmmIBD models IBD along the chromosome and needs linkage
  structure. ADMIXTURE wants the opposite, so this is an easy mistake to make.
  `clonality.ibd_ld_prune` is asserted `false` at parse time and the workflow
  refuses to build if it is ever set true.
- **More markers, but not indiscriminately.** Very low MAF and high-missingness
  markers add noise rather than information.

Non-nuclear contigs (`PKNH_MIT_v2`, `new_API_strain_A1_H.1`) are excluded:
organellar genomes are uniparentally inherited and break hmmIBD's
along-chromosome segment model.

### L2 — Population structure for allele frequencies
**Import cluster *labels* from `data/metadata/admix_clusters.tsv`.** hmmIBD uses
population allele frequencies, so IBD must be computed within genetically
homogeneous groups. Pooling Indonesia + Sabah into one population inflates
apparent IBD between samples from the same region — population structure
masquerading as relatedness.

Importing a cluster *label* is a far lower-stakes import than importing a
pairwise IBD *result*: it is a stable per-sample annotation, and it keeps the
regression comparison apples-to-apples.

Samples with no label form a **residual group**. It is reported in
`cluster_summary.tsv` but **not clonal-called** — a residual group is not a
homogeneous population, so running hmmIBD on it would manufacture exactly the
artefact L2 exists to prevent. Override with
`clonality.residual_clonal_calls: true` if you want to look anyway.

### L3 — Thinning policy
**`declonalize: false` by default; run both.** One representative per clonal
group, chosen by lowest missingness — the same rule `find_duplicates.R` already
uses for replicates, so there is one tie-break convention, not two. Default off
so the headline v1 numbers stay comparable to what is already documented; the
flag makes the sensitivity analysis one parameter away. When false, rule 6 is a
pass-through copy, so the DAG shape is identical either way.

### L4 — Duplicates vs clones
**Duplicates are removed unconditionally; clones are thinned behind L3.**
Removing exact replicates is data hygiene, not a scientific judgement.

The honest limitation: `find_duplicates.R` matches on **ID pattern**, so it
catches technical replicates with lane suffixes and nothing else. A same-patient
re-sample under a different study code is invisible to it and will surface only
as IBD ≈ 1 — indistinguishable from two patients carrying clonal parasites
without metadata. For *discovery* the distinction does not matter (both should be
thinned for frequency estimation). For the *transmission* interpretation later,
it matters a great deal. Recorded as deviation A2b; not solved here.

---

## Rules

| # | Rule | Does |
|---|---|---|
| — | `build_contig_map` | contig → integer map from the reference `.fai`, nuclear only |
| — | `dedupe_missingness` | PLINK2-shaped `.fam`/`.smiss` for the vendored R script |
| 1 | `remove_duplicates` | `find_duplicates.R` → `duplicates.remove` + deduped cohort |
| 2 | `ibd_variant_set` | L1 variant set from the **full** VCF, with assertions |
| 3 | `ibd_cluster_membership` | keep-list per cluster (checkpoint — cluster list is dynamic) |
| 4 | `run_hmmibd` | per cluster, one thread each |
| 5 | `clonal_groups` | `clonal_clusters.R` over every `hmm_fract.txt` |
| 6 | `declonalised_cohort` | apply L3 → `discovery_samples.declonal.txt` |
| 7 | `fws_cross_check` | `run_moimix.R` + `compare_fws_sources.py` — the A1 verdict |

Stage 03's `build_discovery_snp_vcf` consumes
`outputs/qc/discovery_samples.declonal.txt`. Because Stage 03 computes MAF and
F_MISSING *on the sample list it is given*, the SNP set recomputes automatically
— Snakemake gets the "re-run what requires no clones" behaviour from the
dependency edge for free. No orchestration needed.

### Deviations from `agnostic`'s implementation

Both are environment-driven, not design choices, and neither changes a result:

- **No plink2.** `agnostic` builds the variant set and per-cluster subsets with
  `plink2`; this repo has only plink 1.9, and Stages 01–03 are pure bcftools.
  The variant set is built with bcftools instead, keeping the stage in the same
  idiom as the rest of the repo.
- **Missingness source for the replicate tie-break.** `agnostic` uses plink2
  `.smiss` on its `cleaned` bfile; we synthesise a `.smiss` from Stage 02's
  `callable_fraction.tsv` (`F_MISS = 1 − callable_fraction`). The absolute
  numbers differ; only the *ordering within a replicate pair* is used, so the
  chosen representative is unaffected. This keeps `find_duplicates.R`
  byte-identical and avoids a second pass over a 17 GiB VCF.

---

## Running it twice

Same code, one parameter apart. Outputs are keyed by the parameter
(`outputs/ibd/geno10` vs `outputs/ibd/geno20`) so both coexist and can be diffed.

```bash
# Run A — regression. Matches agnostic exactly; proves the port reproduces
# known-good output.
snakemake -s workflow/Snakefile -c8 outputs/qc/clonal_clusters.tsv \
    --config ibd_max_variant_missing=0.10

# Run B — production, per L1. This is the one that feeds Stage 03.
snakemake -s workflow/Snakefile -c8 outputs/qc/clonal_clusters.tsv \
    --config ibd_max_variant_missing=0.20

# Sensitivity arm — thinning on.
snakemake -s workflow/Snakefile -c8 outputs/qc/discovery_samples.declonal.txt \
    --config declonalize=true
```

If A and B give the same clonal groups, the marker set demonstrably does not
matter here and it has been proven rather than assumed. If they diverge, expect
it in Peninsular/Aceh first — the thinnest cluster.

**Result (2026-08-28): they are byte-identical.** Run B carries 68% more variants
overall (553,825 vs 329,805) and 1.6–1.9× more per cluster, and produces the same
21 pairs, 15 groups and 33 samples, member-for-member. The marker set does not
matter here, and that is now proven.

`outputs/qc/clonal_clusters.tsv` always reflects whichever config ran last;
`outputs/ibd/<tag>/clonal_clusters.tsv` is the run-tagged copy, so the two runs
can be diffed directly:

```bash
diff outputs/ibd/geno10/clonal_clusters.tsv outputs/ibd/geno20/clonal_clusters.tsv
```

### Runtime

hmmIBD is single-threaded and scales with pairs × markers; parallelism is
per-cluster, not per-run. On the 593-sample cohort at 330k markers (Run A):
Peninsular (595 pairs) ~16 s, Mn (5,671 pairs) ~4 min, **Mf (81,810 pairs)
~90 min**. Run B roughly doubles Mf. Expect Mf to dominate wall-clock.

---

## Exit criteria — met 2026-08-28

Run A must reproduce this, measured from `agnostic/outputs/ibd/` (2026-08-17), on
the overlapping samples:

| | expected | Run A result |
|---|---|---|
| Clonal pairs at IBD ≥ 0.95 | 21 (Mf 5, Mn 9, Peninsular 7) | **21 (5, 9, 7)** |
| Clonal groups | 15 — 12 pairs + 3 triples | **15 (12 + 3)** |
| Samples in a clonal group | 33 (Mn 15, Mf 10, Peninsular 8) | **33 (15, 10, 8)** |
| Aceh focal pairs | 7, IBD 0.997–1.000 | **7, IBD 0.988–1.000** |

Pair sets are **identical member-for-member** in all three clusters. Note this
reproduction is on **329,805 markers vs agnostic's 24,641**, built with bcftools
rather than plink2 — so the clonal calls are robust to a 13× marker-density change
*and* a different toolchain, which is a stronger result than a like-for-like rerun
would have given.

Cohort: **602 → 593** (9 replicates) → **593** (thinning off by default).
Stage 03 on the 593 gives **51,251 SNPs, up from 49,803** — the upstream MAF-gate
effect predicted in deviation A2, now measured.

Cross-check individual pairs against
`agnostic/outputs/ibd/focal_Peninsular_clones.tsv` — the two Aceh triples there are the
most-scrutinised result in the whole cohort.

Also required:
- `fract_sites_IBD` histogram: bulk near 0, small spike near 1, no bulk in the
  middle (`reports/figures/ibd_fract_distribution_*.png`).
- Duplicate count reported, dropped IDs in the log.
- A1 verdict reported as a number.
- Both cohort sizes stated: post-dedupe and post-thinning.

## Stop conditions

These are enforced in code, not left to a reader's judgement:

| Condition | Enforced where |
|---|---|
| Run A misses 21 pairs by more than ~2 → the port is wrong. **Do not tune thresholds to hit the target.** | manual check against this README |
| Any sample crosses the Fws 0.95 threshold differently between sources → A1 is real; it invalidates the 602 | `compare_fws_sources.py --fail-on-crossing` |
| Duplicate removal drops > ~5% of the cohort → ID pattern over-matching | `remove_duplicates` shell guard |
| `fract_sites_IBD` has substantial mass in the middle → residual polyclonality or an L2 artefact | `plot_ibd.py` prints mid-range % per cluster |
| Thinning removes > ~20% of the cohort → Sabah Pk is more clonal than the discovery design assumes | `apply_declonalisation.py` exits non-zero |

## Gotchas

- **ID normalisation is not the whole story.** The handoff for this stage
  predicted that exact-ID overlap between the two cohorts (546/602) was "purely
  because of Illumina lane suffixes" and that applying `duplicate_id_pattern`
  before the join would fix it. **It does not.** Normalising moves the overlap
  from 546 to 547 — one sample. The remaining ~46 are genuinely absent from
  `admix_clusters.tsv` and land in the residual group. Both joins are still
  applied (exact first, then normalised), but do not expect the second to
  recover much.
- **The right `.fai` is not the configured one.** `reference.fasta` points at the
  PlasmoDB-67 assembly, whose contigs are named with `LT*` accessions. The VCF
  uses `ordered_PKNH_*_v2`. The contig map is built from
  `data/reference/strain_preDB_version/strain_A1_H.1.Icor.fasta.fai`, which is
  the one that matches. Using the configured `.fai` yields an empty map and a
  cryptic failure inside `genotype_table.R`.
- **No LD pruning.** Stated twice deliberately; asserted in code.
- **Ploidy.** Our VCF is diploid-called (`0/0`, `1/1`) with a 0.48% het rate and
  6.5% per-call missingness — fine for hmmIBD, which expects exactly this shape.
- **`hmmIBD` is built from source** into the project env by
  `envs/hmmIBD-src/hmmIBD.c`; it is not in the conda lockfile.
- **Do not modify anything under `agnostic/`.** Copy out only.
