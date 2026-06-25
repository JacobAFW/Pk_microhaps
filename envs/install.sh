#!/usr/bin/env bash
# install.sh — bootstrap the P. knowlesi microhap panel-discovery environment
#
# One-shot install for collaborators. Installs vvg-box into the project root,
# then populates it with every tool the v1 pipeline needs:
#   - Python: scikit-allel, numpy, pandas, matplotlib, jupyter (Siegel discovery + selection)
#   - R: paneljudge, moimix, ape, igraph, tidyverse (validation + plotting)
#   - External: hmmIBD v3, THE REAL McCOIL v2 (validation, deferred to Stage 06)
#   - CLI: bcftools, htslib, samtools, plink, snakemake, quarto
#
# Idempotent: rerunning after a partial install picks up where it left off.
#
# Usage:
#   cd /path/to/Microhaplotype
#   bash envs/install.sh
#
# VCF handling — the discovery cohort VCF (~17 GiB) is the project's biggest
# dependency. Three modes, in order of preference:
#   1. POPGEN_PIPELINE_PATH set → symlink data/vcf/ to that path's data/vcf/
#   2. /path/to/Pop-gen_pipeline/ found at the default sibling path → symlink
#   3. Neither → leave data/vcf/ empty and print copy-mode instructions
#
# Scope: all installs stay inside the project. Nothing is written to /usr/local,
# /opt, ~/.zshrc, or any system directory.

set -euo pipefail

#---------------------------------------------------------------------------
# Re-exec under bash >= 4 if needed (vvg-box source uses [[ -v ... ]],
# associative arrays, etc. — none of which work in macOS-stock bash 3.2).
#---------------------------------------------------------------------------
if (( ${BASH_VERSINFO[0]:-0} < 4 )); then
  for newer_bash in /opt/homebrew/bin/bash /usr/local/bin/bash; do
    if [[ -x "$newer_bash" ]]; then
      echo "==> Re-executing under $newer_bash (current bash ${BASH_VERSION} too old)"
      exec "$newer_bash" "$0" "$@"
    fi
  done
  echo "ERROR: bash >= 4 required (have ${BASH_VERSION}). Install with: brew install bash" >&2
  exit 1
fi

#---------------------------------------------------------------------------
# Locate the project root
#---------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

if [[ ! -f "$PROJECT_ROOT/INSTRUCTIONS.md" ]]; then
  echo "ERROR: could not locate project root (no INSTRUCTIONS.md at $PROJECT_ROOT)"
  exit 1
fi

cd "$PROJECT_ROOT"
echo "==> Project root: $PROJECT_ROOT"
echo "==> Architecture: $(uname -s)-$(uname -m)"

#---------------------------------------------------------------------------
# macOS coreutils shim (vvg-box expects GNU ln/readlink semantics)
#---------------------------------------------------------------------------
if [[ "$(uname -s)" == "Darwin" ]]; then
  for gnubin in /opt/homebrew/opt/coreutils/libexec/gnubin /usr/local/opt/coreutils/libexec/gnubin; do
    if [[ -d "$gnubin" ]]; then
      export PATH="$gnubin:$PATH"
      echo "==> Prepending GNU coreutils to PATH: $gnubin"
      break
    fi
  done
  # Put Homebrew bin before /bin so subprocess `bash` calls (e.g. the curl-bash
  # for vvg-box's installer) pick up Homebrew's bash 5 rather than /bin/bash 3.2.
  for hb in /opt/homebrew/bin /usr/local/bin; do
    if [[ -x "$hb/bash" ]]; then
      export PATH="$hb:$PATH"
      echo "==> Prepending $hb to PATH (bash $($hb/bash --version | head -1 | awk '{print $4}'))"
      break
    fi
  done
  if ! ln --help 2>&1 | grep -q -- '--relative'; then
    echo "ERROR: macOS detected but GNU coreutils not found."
    echo "  Install with: brew install coreutils"
    exit 1
  fi
fi

#---------------------------------------------------------------------------
# Step 1: install vvg-box if not already present
#
# Keep the install fully project-scoped: override MAMBA_ROOT_PREFIX so the
# bundled micromamba env lands inside the project rather than the user's
# personal ~/mamba root. Without this, vvg-box's installer collides with
# any pre-existing user-level mamba root that happens to have a `vvg-box`
# env name registered.
#
# Important: the path must NOT be ./vvg-box — vvg-box's installer also
# git-clones its own source tree into ./vvg-box, so making that the mamba
# root collides with the env path ./vvg-box/envs/vvg-box.
#---------------------------------------------------------------------------
export MAMBA_ROOT_PREFIX="$PROJECT_ROOT/.vvg-mamba"

# Treat install as "already done" only if the mamba env activate script exists.
# Skeleton dirs alone (created by a partial/aborted run) are not enough.
VVGBOX_ENV_ACTIVATE="$MAMBA_ROOT_PREFIX/envs/vvg-box/bin/activate"
if [[ ! -f "$VVGBOX_ENV_ACTIVATE" \
   && ! -f "$PROJECT_ROOT/activate" ]]; then
  echo "==> Installing vvg-box into $PROJECT_ROOT (MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX)"
  # If a partial install left behind a half-skeleton or empty env, clear it so
  # the vvg-box installer can start fresh and not trip its own overwrite prompt.
  if [[ -d "$PROJECT_ROOT/vvg-box" ]]; then
    echo "==> Removing partial vvg-box skeleton from previous run"
    rm -rf "$PROJECT_ROOT/vvg-box"
  fi
  if [[ -d "$MAMBA_ROOT_PREFIX" ]]; then
    echo "==> Removing partial mamba root from previous run"
    rm -rf "$MAMBA_ROOT_PREFIX"
  fi
  # vvg-box's installer is a bash-only script (uses repeat() etc. that zsh refuses).
  bash <(curl -L https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh)
else
  echo "==> vvg-box already installed — skipping bootstrap"
fi

#---------------------------------------------------------------------------
# Step 2: activate vvg-box for this install session
#---------------------------------------------------------------------------
if [[ -f "$PROJECT_ROOT/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/activate"
elif [[ -f "$PROJECT_ROOT/vvg-box/bin/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/vvg-box/bin/activate"
elif [[ -f "$PROJECT_ROOT/bin/activate" ]]; then
  ACTIVATE="$PROJECT_ROOT/bin/activate"
else
  echo "ERROR: could not find vvg-box activation script. Tried:"
  echo "  $PROJECT_ROOT/activate"
  echo "  $PROJECT_ROOT/vvg-box/bin/activate"
  echo "  $PROJECT_ROOT/bin/activate"
  exit 1
fi

# vvg-box's activate + bashrc.d trip set -u (unbound vars) and history -r.
set +eu
# shellcheck disable=SC1090
source "$ACTIVATE"
set -eu
echo "==> Activated vvg-box from $ACTIVATE"

# Record activation path for future sessions
cat > "$PROJECT_ROOT/envs/activate.sh" <<EOF
#!/usr/bin/env bash
# Source this file to activate the pipeline environment.
# Generated by envs/install.sh on $(date).
source "$ACTIVATE"
EOF
chmod +x "$PROJECT_ROOT/envs/activate.sh"

#---------------------------------------------------------------------------
# Step 3: pick the conda-family installer vvg-box exposed
#---------------------------------------------------------------------------
if command -v mamba &>/dev/null; then
  INSTALLER="mamba install -y -c bioconda -c conda-forge"
elif command -v micromamba &>/dev/null; then
  INSTALLER="micromamba install -y -c bioconda -c conda-forge"
elif command -v conda &>/dev/null; then
  INSTALLER="conda install -y -c bioconda -c conda-forge"
else
  echo "ERROR: no mamba/micromamba/conda on PATH after vvg-box activation"
  exit 1
fi
echo "==> Installer: $INSTALLER"

#---------------------------------------------------------------------------
# Step 4: CLI tools (match Pop-gen Stage-1/2 pins where they apply)
#---------------------------------------------------------------------------
echo "==> Installing CLI tools (bcftools, htslib, samtools, plink, bedtools, snakemake, quarto)"
$INSTALLER \
  bcftools=1.21 \
  htslib=1.21 \
  samtools \
  plink \
  bedtools \
  snakemake-minimal \
  quarto

#---------------------------------------------------------------------------
# Step 5: Python stack for the Siegel discovery framework port
#   Pins follow vivax-mhaps/pyproject.toml: numpy<=1.23.5, pandas^2.2,
#   scikit-allel<=1.3.7, matplotlib^3.8.2. malariagen_data is NOT installed —
#   we feed our own Pk VCF rather than pulling from MalariaGEN cloud.
#---------------------------------------------------------------------------
echo "==> Installing Python stack (Siegel discovery deps)"
$INSTALLER \
  python=3.11 \
  "numpy<=1.23.5" \
  "pandas>=2.2,<3.0" \
  "scikit-allel<=1.3.7" \
  "matplotlib>=3.8,<4.0" \
  scikit-learn=1.2.0 \
  jupyter \
  pyvcf3 \
  pysam

#---------------------------------------------------------------------------
# Step 6: R stack for validation + plotting (Stage 06)
#---------------------------------------------------------------------------
echo "==> Installing R stack via mamba (CRAN/bioconda where available)"
$INSTALLER \
  r-base=4.3.1 \
  r-tidyverse \
  r-ape \
  r-igraph \
  r-data.table \
  r-biocmanager \
  r-remotes

# paneljudge — not on CRAN, GitHub install (aimeertaylor/paneljudge per Siegel methods)
echo "==> Installing paneljudge from GitHub (aimeertaylor/paneljudge)"
Rscript -e '
  if (!requireNamespace("paneljudge", quietly = TRUE)) {
    remotes::install_github("aimeertaylor/paneljudge", upgrade = "never")
  } else {
    message("paneljudge already installed — skipping")
  }
'

# moimix — needed if we want to recompute Fws standalone (we will reuse Pop-gen's
# output by default, but having moimix here means the project can run without Pop-gen)
echo "==> Installing moimix from GitHub (bahlolab/moimix)"
Rscript -e '
  if (!requireNamespace("moimix", quietly = TRUE)) {
    remotes::install_github("bahlolab/moimix", upgrade = "never")
  } else {
    message("moimix already installed — skipping")
  }
'

#---------------------------------------------------------------------------
# Step 7: hmmIBD + THE REAL McCOIL — DEFERRED to Stage 06 install pass
#---------------------------------------------------------------------------
# Both are needed only for validation (Stage 06). Building from source adds
# ~5 min and isn't needed for Stages 01–05. Re-run install.sh with
# INSTALL_VALIDATION=1 to add them.
if [[ "${INSTALL_VALIDATION:-0}" == "1" ]]; then
  echo "==> Building hmmIBD v3 from source"
  HMMIBD_SRC="$PROJECT_ROOT/external/clones/hmmIBD"
  if [[ ! -d "$HMMIBD_SRC" ]]; then
    git clone --depth 1 https://github.com/glipsnort/hmmIBD.git "$HMMIBD_SRC"
  fi
  ( cd "$HMMIBD_SRC" && make )
  ln -sf "$HMMIBD_SRC/hmmIBD" "$PROJECT_ROOT/envs/bin/hmmIBD" 2>/dev/null || true

  echo "==> Building THE REAL McCOIL v2 from source"
  RMC_SRC="$PROJECT_ROOT/external/clones/THEREALMcCOIL"
  if [[ ! -d "$RMC_SRC" ]]; then
    git clone --depth 1 https://github.com/EPPIcenter/THEREALMcCOIL.git "$RMC_SRC"
  fi
  # build steps depend on the variant (categorical vs proportional); see
  # CONTEXT-software.md for the chosen build.
  echo "==> NOTE: REAL McCOIL build is variant-specific — see CONTEXT-software.md"
else
  echo "==> Skipping validation tools (hmmIBD, REAL McCOIL). Re-run with INSTALL_VALIDATION=1 when ready for Stage 06."
fi

#---------------------------------------------------------------------------
# Step 8: VCF + reference data — symlink from Pop-gen, or instruct copy
#---------------------------------------------------------------------------
echo
echo "==> Setting up data/vcf/ and data/reference/"

POPGEN_DEFAULT="$( cd "$PROJECT_ROOT/.." && pwd )/Pop-gen_pipeline"
POPGEN_PATH="${POPGEN_PIPELINE_PATH:-$POPGEN_DEFAULT}"

if [[ -d "$POPGEN_PATH/data/vcf" && -f "$POPGEN_PATH/data/vcf/merged_popgen.vcf.gz" ]]; then
  echo "==> Pop-gen_pipeline found at: $POPGEN_PATH"
  echo "==> Symlinking data/vcf/ and data/reference/ from Pop-gen"

  # VCF
  ln -sfn "$POPGEN_PATH/data/vcf/merged_popgen.vcf.gz" \
          "$PROJECT_ROOT/data/vcf/merged_popgen.vcf.gz"
  if [[ -f "$POPGEN_PATH/data/vcf/merged_popgen.vcf.gz.csi" ]]; then
    ln -sfn "$POPGEN_PATH/data/vcf/merged_popgen.vcf.gz.csi" \
            "$PROJECT_ROOT/data/vcf/merged_popgen.vcf.gz.csi"
  fi

  # Reference (PlasmoDB-67 fasta + indices, PlasmoDB-68 GFF)
  for f in PlasmoDB-67_PknowlesiA1H1_Genome.fasta \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.fai \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.amb \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.ann \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.bwt \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.pac \
           PlasmoDB-67_PknowlesiA1H1_Genome.fasta.sa; do
    src="$POPGEN_PATH/data/reference/PlasmoDB_version/$f"
    [[ -f "$src" ]] && ln -sfn "$src" "$PROJECT_ROOT/data/reference/$f"
  done
  for f in PlasmoDB-68_PknowlesiA1H1.gff PlasmoDB-68_PknowlesiA1H1_GO.gaf; do
    src="$POPGEN_PATH/data/reference/PlasmoDB_version/$f"
    [[ -f "$src" ]] && ln -sfn "$src" "$PROJECT_ROOT/data/reference/$f"
  done
  if [[ -f "$POPGEN_PATH/data/reference/regions_to_mask.list" ]]; then
    ln -sfn "$POPGEN_PATH/data/reference/regions_to_mask.list" \
            "$PROJECT_ROOT/data/reference/regions_to_mask.list"
  fi

  echo "==> Symlinks in place. data/vcf/ and data/reference/ now point at Pop-gen."

  # Persist the resolved POPGEN_PIPELINE_PATH so future shells (and snakemake)
  # find the Fws table without re-supplying the env var on every command.
  if ! grep -q '^export POPGEN_PIPELINE_PATH=' "$PROJECT_ROOT/envs/activate.sh"; then
    echo "export POPGEN_PIPELINE_PATH=\"$POPGEN_PATH\"" >> "$PROJECT_ROOT/envs/activate.sh"
    echo "==> Recorded POPGEN_PIPELINE_PATH in envs/activate.sh"
  fi
else
  cat <<EOF

==> Pop-gen_pipeline NOT found.
    Tried: $POPGEN_PATH
    Either:
      (a) clone or copy Pop-gen_pipeline next to this directory and re-run, or
      (b) set POPGEN_PIPELINE_PATH=/path/to/your/Pop-gen_pipeline and re-run, or
      (c) drop the discovery VCF + reference manually into:
          $PROJECT_ROOT/data/vcf/merged_popgen.vcf.gz (+ .csi)
          $PROJECT_ROOT/data/reference/PlasmoDB-67_PknowlesiA1H1_Genome.fasta (+ .fai, BWA indices)
          $PROJECT_ROOT/data/reference/PlasmoDB-68_PknowlesiA1H1.gff
          $PROJECT_ROOT/data/reference/regions_to_mask.list

EOF
fi

#---------------------------------------------------------------------------
# Step 9: lockfile
#---------------------------------------------------------------------------
LOCKFILE="$PROJECT_ROOT/envs/environment.lock.yaml"
echo "==> Writing lockfile to $LOCKFILE"
if command -v micromamba &>/dev/null; then
  micromamba env export --explicit > "$LOCKFILE" 2>/dev/null || micromamba env export > "$LOCKFILE"
elif command -v mamba &>/dev/null; then
  mamba env export > "$LOCKFILE"
fi

#---------------------------------------------------------------------------
# Step 10: sanity check
#---------------------------------------------------------------------------
echo
echo "==> Sanity-check installed tools:"
for tool in bcftools samtools plink Rscript snakemake quarto python jupyter; do
  if command -v "$tool" &>/dev/null; then
    ver="$("$tool" --version 2>&1 | head -1 || true)"
    printf "    %-12s  %s\n" "$tool" "$ver"
  else
    printf "    %-12s  NOT FOUND\n" "$tool"
  fi
done
echo
echo "==> Python module check:"
python -c "import allel, numpy, pandas, sklearn, matplotlib; \
  print(f'  scikit-allel  {allel.__version__}'); \
  print(f'  numpy         {numpy.__version__}'); \
  print(f'  pandas        {pandas.__version__}'); \
  print(f'  scikit-learn  {sklearn.__version__}'); \
  print(f'  matplotlib    {matplotlib.__version__}')"

echo
echo "==> Done. Activate in future sessions with:"
echo "     source envs/activate.sh"
echo
echo "==> Validate the Snakemake DAG with:"
echo "     snakemake -n"
echo
echo "==> Check CONTEXT-software.md for tool-version deviations from Pop-gen / paper."
