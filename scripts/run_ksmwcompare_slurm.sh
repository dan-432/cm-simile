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
# REPO_DIR is hardcoded (not derived from $SLURM_SUBMIT_DIR/$PWD) because
# that's the exact failure mode that broke previous runs: sbatch can be
# launched from inside scripts/, so $PWD resolves to .../cm-simile/scripts
# rather than .../cm-simile, which leads to a failed Go build from the wrong
# working directory. Keep this path in sync with the repo location.
REPO_DIR="/projects/health_sciences/bms/biochemistry/gardner_group/daniel/cm-simile"
RFAM_CM="${REPO_DIR}/data/Rfam.cm"
CLANIN="${REPO_DIR}/data/Rfam.clanin"
BINARY="${REPO_DIR}/bin/ksmwcompare"
POP_DIR="${REPO_DIR}/results/clan_population_${SLURM_JOB_ID}"
BENCH_DIR="${REPO_DIR}/results/ksmwcompare_benchmark_${SLURM_JOB_ID}"
N_CLANS=8
N_SINGLETON_DISTRACTORS=20
N_CROSSCLAN_DISTRACTORS=20
NPOOL=1000
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

mkdir -p "$POP_DIR" "$BENCH_DIR"

# -------------------------------------------------------------------------
# Step 1: build the clan population (CM file + clan annotation TSV).
# Same construction script as the other benchmark workflows -- point POP_DIR
# at an earlier run's output and comment this block out if you want to reuse
# the same population across multiple benchmark runs.
# -------------------------------------------------------------------------
python3 scripts/constructtestpopulation.py \
  --rfam-cm "$RFAM_CM" \
  --clanin "$CLANIN" \
  --n-clans "$N_CLANS" \
  --n-singleton-distractors "$N_SINGLETON_DISTRACTORS" \
  --n-crossclan-distractors "$N_CROSSCLAN_DISTRACTORS" \
  --seed "$SEED" \
  --out "$POP_DIR"

# -------------------------------------------------------------------------
# Step 2: run the ksmwcompare benchmark over that population.
# -------------------------------------------------------------------------
python3 scripts/ksmwcomparebenchmark.py \
  --cm-file "$POP_DIR/population.cm" \
  --binary "$BINARY" \
  --npool "$NPOOL" \
  --n-null-pools "$N_NULL_POOLS" \
  --seed "$SEED" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --out "$BENCH_DIR/raw_results.tsv"

echo "[slurm] done, $(date). Population in $POP_DIR, results in $BENCH_DIR"
