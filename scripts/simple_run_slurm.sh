#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cmsimile_ksmw
#SBATCH --output=cmsimile_%j.log
#SBATCH --time=35:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate cm-simile

set -euo pipefail
# -------------------------------------------------------------------------
# Config -- edit these per run
# -------------------------------------------------------------------------

REPO_DIR="/projects/health_sciences/bms/biochemistry/gardner_group/daniel/cm-simile"
CM="${REPO_DIR}/data/all_iteratively_expanded_rf00050_dediscovered_and_references.cm"
BINARY="${REPO_DIR}/bin/ksmwcompare"
BENCH_DIR="${REPO_DIR}/results/ksmwcompare_diagnostic_${SLURM_JOB_ID}"
NPOOL=300
N_NULL_POOLS=1
SEED=42
# -------------------------------------------------------------------------
# Environment
# -------------------------------------------------------------------------

cd "$REPO_DIR"

# Fail fast with a clear message rather than a cryptic Go toolchain error
# if REPO_DIR is ever wrong again (repo moved, script copied elsewhere,
# submitted from an unexpected location, etc.).
if [ ! -f "$REPO_DIR/go.mod" ]; then
    echo "[slurm] ERROR: no go.mod found at $REPO_DIR -- REPO_DIR is wrong. Check the hardcoded path at the top of this script." >&2
    exit 1
fi

echo "[slurm] job $SLURM_JOB_ID starting on $(hostname), $SLURM_CPUS_PER_TASK cpus, $(date)"
# -------------------------------------------------------------------------
# Build the binary if it isn't already there. Cheap and idempotent, so
# it's safe to leave in even once bin/ is committed/cached between runs.
# -------------------------------------------------------------------------
if [ ! -x "$BINARY" ]; then
    echo "[slurm] building $BINARY ..."
    mkdir -p "$(dirname "$BINARY")"
    go build -o "$BINARY" src/ksmwcompare.go
fi


# -------------------------------------------------------------------------
# run ksmwcompare
# -------------------------------------------------------------------------

mkdir -p "$BENCH_DIR"

"$BINARY" -cm-file "$CM" \
    -out "$BENCH_DIR/ksmwcompare-diagnostic-analysis.tsv" \
    --npool "$NPOOL" \
    --n-null-pools "$N_NULL_POOLS" \
    -seed "$SEED" \
    -threads "$SLURM_CPUS_PER_TASK" \
    -query "FMN_on_micro"


echo "[slurm] done, $(date). Results in $BENCH_DIR"
