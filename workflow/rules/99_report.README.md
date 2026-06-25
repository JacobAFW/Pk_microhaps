# 99_report — Quarto book

**Status:** implemented (HANDOFF #03, 2026-05-12). Grows by one chapter per completed stage. Not yet rendered.

## What this stage does
Renders the project Quarto **book** from per-stage chapters in `reports/qmd/` to a multi-page HTML site at `reports/_book/`. The book is re-rendered after each stage completes — readers always see a current view through whichever stages are done.

Book config: `_quarto.yml` (project root). Execute-dir is set to project root so chapter code cells read outputs with project-root-relative paths.

## Chapters (in book order)
| Chapter | Status |
|---|---|
| `reports/qmd/index.qmd` | Implemented. Project overview + reading guide. |
| `reports/qmd/01_data_prep.qmd` | **Stub** — back-fill pending before HANDOFF #04. Has key-numbers block + summary. |
| `reports/qmd/02_sample_qc.qmd` | Implemented. |
| `reports/qmd/03_snp_qc.qmd` | Implemented. |
| `reports/qmd/04_window_scan.qmd` | Not yet (stage not yet ported). |
| `reports/qmd/05_panel_select.qmd` | Not yet. |
| `reports/qmd/06_validate.qmd` | Not yet. |

## Inputs (current)
- `_quarto.yml`
- `reports/qmd/{index,01_data_prep,02_sample_qc,03_snp_qc}.qmd`
- All Stage 01–03 outputs the chapters embed (BED, Fws table, sample_metadata, callable_fraction, discovery list, snps.discovery VCF + summary, all five figures).

## Outputs
- `reports/_book/index.html` — book entry point.
- Sibling chapter HTMLs and supporting assets under `reports/_book/`.

## Rules
- `render_report` — runs `quarto render` from the project root; depends on every chapter source + every output the chapters embed. Snakemake re-renders when any input changes.

## Audience
Two tiers:
- **Jacob** (technical): the source of truth for what each Claude Code run actually produced, the parameters that were used, and what the figures look like.
- **The boss / collaborators** (boss-readable): can skim the per-chapter summary blocks and the embedded outputs without expanding code cells.

Code cells are folded by default (`code-fold: true`) — click to expand.

## Adding a chapter
When a new stage lands:
1. Add `reports/qmd/NN_stage.qmd` matching the existing chapter shape (Summary → Key numbers → Methods → Parameters → Outputs → Figures → Decisions → Caveats → What's next).
2. Add it to `_quarto.yml` under `book.chapters`.
3. Add its source + the stage's outputs to `99_report.smk` `render_report.input:`.
4. Re-render via `snakemake reports/_book/index.html`.

## Audit / reproducibility
- `execute.freeze: auto` in `_quarto.yml` caches chapter renders so an unchanged chapter doesn't re-execute. To force a clean re-render, delete `reports/qmd/_freeze/` or pass `--freeze-mode auto-redo` to the render command.
- The render rule's log lives at `logs/99_report/render_report.log` with the Quarto stdout / stderr.
