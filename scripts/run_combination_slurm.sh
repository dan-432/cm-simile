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
CM="${REPO_DIR}/data/all_clusters_an_mafft_aligned_combined_and_true_model.cm"
BINARY="${REPO_DIR}/bin/ksmwcompare"
BENCH_DIR="${REPO_DIR}/results/ksmwcompare_combine_${SLURM_JOB_ID}"
NPOOL=300
N_NULL_POOLS=1
COMBINE_THRESHOLD_KS=0.8
COMBINE_TOP_PROPORTION=0.7
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

# mkdir -p "$BENCH_DIR/combined_models"


# -------------------------------------------------------------------------
# run the ksmwcompare with combination mode over that population.
# -------------------------------------------------------------------------

# "$BINARY" -cm-file "$CM" \
#     -npool "$NPOOL" -combine -combine-threshold-ks "$COMBINE_THRESHOLD_KS" -combine-top-proportion "$COMBINE_TOP_PROPORTION" \
#     -out "$BENCH_DIR/ksmw-combine.tsv" \
#     -combine-outdir "$BENCH_DIR/combined_models" \
#     -n-null-pools "$N_NULL_POOLS" \
#     -seed "$SEED" \
#     -threads "$SLURM_CPUS_PER_TASK"

"$BINARY" -cm-file "$CM" \
    -npool "$NPOOL" \
    -out "$BENCH_DIR/ksmw-combined-analysis.tsv" \
    -n-null-pools "$N_NULL_POOLS" \
    -seed "$SEED" \
    -threads "$SLURM_CPUS_PER_TASK"


echo "[slurm] done, $(date). Results in $BENCH_DIR"
