#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cmsimile_dan
#SBATCH --output=cmsimile_%j.log
#SBATCH --time=35:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32GB                 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate cm-simile


set -euo pipefail
# -------------------------------------------------------------------------
# Config -- edit these per run
# -------------------------------------------------------------------------
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
RFAM_CM="${REPO_DIR}/data/Rfam.cm"
CLANIN="${REPO_DIR}/data/Rfam.clanin"
BINARY="${REPO_DIR}/bin/minmaxcmcompare"
POP_DIR="${REPO_DIR}/results/clan_population_${SLURM_JOB_ID}"
BENCH_DIR="${REPO_DIR}/results/minmax_benchmark_${SLURM_JOB_ID}"
N_CLANS=8
N_SINGLETON_DISTRACTORS=20
N_CROSSCLAN_DISTRACTORS=20
NGENERATIONS=1
NPOPULATION=200
SEED=42
# -------------------------------------------------------------------------
# Environment
# -------------------------------------------------------------------------

cd "$REPO_DIR"
echo "[slurm] job $SLURM_JOB_ID starting on $(hostname), $SLURM_CPUS_PER_TASK cpus, $(date)"
# -------------------------------------------------------------------------
# Build the binary if it isn't already there. Cheap and idempotent, so
# it's safe to leave in even once bin/ is committed/cached between runs --
# saves a "forgot to rebuild after a code change" class of mistake.
#
# NB: assuming the source file is src/minmaxcmcompare.go to match the
# binary name -- adjust if the Go source is named differently.
# -------------------------------------------------------------------------
if [ ! -x "$BINARY" ]; then
echo "[slurm] building $BINARY ..."
mkdir -p "$(dirname "$BINARY")"
go build -o "$BINARY" src/minmaxcmcompare.go
fi
mkdir -p "$POP_DIR" "$BENCH_DIR"
# -------------------------------------------------------------------------
# Step 1: build the clan population (CM file + clan-annotation TSV).
# Deliberately kept separate from the benchmark run below, so the same
# population can be reused across multiple benchmark/parameter sweeps
# without rebuilding it each time -- for a sweep, comment this block out
# on subsequent submissions and point POP_DIR at the earlier run's dir.
# -------------------------------------------------------------------------
python3 constructtestpopulation.py \
--rfam-cm "$RFAM_CM" \
--clanin "$CLANIN" \
--n-clans "$N_CLANS" \
--n-singleton-distractors "$N_SINGLETON_DISTRACTORS" \
--n-crossclan-distractors "$N_CROSSCLAN_DISTRACTORS" \
--seed "$SEED" \
--out "$POP_DIR"
# -------------------------------------------------------------------------
# Step 2: run the minmaxcompare benchmark over that population.
# -------------------------------------------------------------------------
python3 minmaxbenchmark.py \
--cm-file "$POP_DIR/population.cm" \
--binary "$BINARY" \
--generations "$NGENERATIONS" \
--pop-size "$NPOPULATION" \
--seed "$SEED" \
--threads "$SLURM_CPUS_PER_TASK" \
--out "$BENCH_DIR/raw_results.tsv"
echo "[slurm] done, $(date). Population in $POP_DIR, results in $BENCH_DIR"
