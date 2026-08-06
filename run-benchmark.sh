#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=cmsimile_dan
#SBATCH --output=cmsimile_%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
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
BINARY="${REPO_DIR}/bin/bitscoresumdifference"
OUTDIR="${REPO_DIR}/results/clan_benchmark_${SLURM_JOB_ID}"
 
N_CLANS=8
N_SINGLETON_DISTRACTORS=20
N_CROSSCLAN_DISTRACTORS=20
NPOOL=100
NSHUFFLES=200
SEED=42
 
# -------------------------------------------------------------------------
# Environment
# -------------------------------------------------------------------------
source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate cmsimile
 
cd "$REPO_DIR"
 
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
 
mkdir -p "$OUTDIR"
 
# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------
python3 scripts/run_clan_benchmark.py \
    --rfam-cm "$RFAM_CM" \
    --clanin "$CLANIN" \
    --binary "$BINARY" \
    --n-clans "$N_CLANS" \
    --n-singleton-distractors "$N_SINGLETON_DISTRACTORS" \
    --n-crossclan-distractors "$N_CROSSCLAN_DISTRACTORS" \
    --npool "$NPOOL" \
    --nshuffles "$NSHUFFLES" \
    --seed "$SEED" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --out "$OUTDIR"
 
echo "[slurm] done, $(date). Results in $OUTDIR"