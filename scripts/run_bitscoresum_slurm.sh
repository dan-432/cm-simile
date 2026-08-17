#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cmsimile_dan
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
# that's exactly what broke last run: sbatch was invoked from inside
# scripts/, so ${SLURM_SUBMIT_DIR:-$PWD} silently resolved to
# .../cm-simile/scripts instead of .../cm-simile, `cd`'d there, and `go
# build` -- finding no go.mod in that directory -- fell back to legacy
# GOPATH-style import resolution, producing the "package ... is not in
# std" error (which is really Go saying "I couldn't find a module here
# and guessed you meant a stdlib import path", not a real problem with
# bitscoresumdifference.go itself). Update this path if the repo moves.
REPO_DIR="/projects/health_sciences/bms/biochemistry/gardner_group/daniel/cm-simile"
RFAM_CM="${REPO_DIR}/data/Rfam.cm"
CLANIN="${REPO_DIR}/data/Rfam.clanin"
BINARY="${REPO_DIR}/bin/bitscoresumdifference"
POP_DIR="${REPO_DIR}/results/clan_population_${SLURM_JOB_ID}"
BENCH_DIR="${REPO_DIR}/results/bitscoresum_benchmark_${SLURM_JOB_ID}"
N_CLANS=8
N_SINGLETON_DISTRACTORS=20
N_CROSSCLAN_DISTRACTORS=20
NPOOL=100
NSHUFFLES=200
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
# it's safe to leave in even once bin/ is committed/cached between runs --
# saves a "forgot to rebuild after a code change" class of mistake.
# -------------------------------------------------------------------------
if [ ! -x "$BINARY" ]; then
echo "[slurm] building $BINARY ..."
mkdir -p "$(dirname "$BINARY")"
go build -o "$BINARY" src/bitscoresumdifference.go
fi
mkdir -p "$POP_DIR" "$BENCH_DIR"
# -------------------------------------------------------------------------
# Step 1: build the clan population (CM file + clan-annotation TSV).
# Same construction script as the minmaxcmcompare workflow -- point
# POP_DIR at an earlier run's output (and comment this block out) to
# benchmark both algorithms against the IDENTICAL population rather than
# two independently-sampled ones, which matters if you want a fair
# side-by-side comparison between the two algorithms later.
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
# Step 2: run the bitscoresumdifference benchmark over that population.
# -------------------------------------------------------------------------
python3 scripts/bitscoresumdifferencebenchmark.py \
--cm-file "$POP_DIR/population.cm" \
--binary "$BINARY" \
--npool "$NPOOL" \
--nshuffles "$NSHUFFLES" \
--seed "$SEED" \
--threads "$SLURM_CPUS_PER_TASK" \
--out "$BENCH_DIR/raw_results.tsv"
echo "[slurm] done, $(date). Population in $POP_DIR, results in $BENCH_DIR"