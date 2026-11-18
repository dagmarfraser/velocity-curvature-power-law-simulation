#!/bin/bash
#SBATCH --job-name=v028h
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/slurm_v028_%j.out
#SBATCH --mail-type=END,FAIL
#
# runV028HPC_v003.sh  (2026-11-16, v003 adds --harmonics)  Batch job for runLoopClosureFftnoise_v028 (supersedes runV024HPC_v001.sh): end-anchored template, detectability rule primary; --rule budget for the 1% arm.
# Submit FROM the RDS _N48 src/ folder (docs/README_BlueBEARBatch_v001.md), always with --account=fraserds-mpo-evaluation:
#   primary, core + out-of-domain arms:
#     sbatch --account=fraserds-mpo-evaluation runV028HPC_v003.sh "Hickman PLAC" "Hickman HALO" "Cook CTRL" "Cook ASD" Zarandi Dhieb
#     sbatch --account=fraserds-mpo-evaluation runV028HPC_v003.sh FraserN48
#   shapeXu sensitivity arm:      ... runV028HPC_v003.sh --model shapedXuSeg <datasets>
#   1% budget arm:                ... runV028HPC_v003.sh --rule budget <datasets>
#   seed-stability rerun:         ... runV028HPC_v003.sh --seed 271828 --reps 50 <datasets>   (own output names, no overwrite)
#   harmonic-count sensitivity:   ... runV028HPC_v003.sh --harmonics fc <datasets>   (fc = harmonics up to the 10 Hz cutoff; output names add _nHfc)
#   acceptance (small grid):      sbatch --account=... --cpus-per-task=16 --mem=64G --time=1:00:00 runV028HPC_v003.sh --acceptance "Cook CTRL"
# Watch: squeue -u $USER; tail -f logs/slurm_v028_<jobid>.out; sacct -j <jobid>; seff <jobid>. Resume: submit the same line again.
set -euo pipefail

module purge; module load bluebear
module load bear-apps/2023a
module load MATLAB/2024b-r6

ACCEPT=false; MODEL=endMatched; SEED=1729; REPS=200; RULE=detect; HARM=4
while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
        --acceptance) ACCEPT=true; shift ;;
        --model) MODEL="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --reps) REPS="$2"; shift 2 ;;
        --rule) RULE="$2"; shift 2 ;;
        --harmonics) HARM="$2"; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
done
if [[ $# -lt 1 ]]; then echo "usage: sbatch runV028HPC_v003.sh [--acceptance] [--model endMatched|shapedXuSeg] [--seed N] [--reps N] [--rule detect|budget] [--harmonics 4|fc|N] <dataset> ..." >&2; exit 2; fi

export TMPDIR="${TMPDIR:-/tmp}/matlab_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
DS=""; for d in "$@"; do DS+="\"$d\","; done; DS="[${DS%,}]"
END=$(squeue -h -j "$SLURM_JOB_ID" -o %e)

echo "Job $SLURM_JOB_ID on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cores, datasets $DS, model $MODEL, rule $RULE, seed $SEED, reps $REPS, ends $END, acceptance $ACCEPT"
if $ACCEPT; then
    matlab -batch "runV028HPC_v003($DS, \"\", Acceptance=true)"
else
    matlab -batch "runV028HPC_v003($DS, \"$END\", NoiseModel=\"$MODEL\", Rule=\"$RULE\", Harmonics=\"$HARM\", RngSeed=$SEED, NReps=$REPS)"
fi
echo "MATLAB finished at $(date)"
