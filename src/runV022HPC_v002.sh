#!/bin/bash
#SBATCH --job-name=v022_endMatched
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/slurm_v022_%j.out
#SBATCH --mail-type=END,FAIL
#
# runV022HPC_v002.sh  (2026-11-15; v002 calls runV022HPC_v002.m, whose acceptance mode uses checkV022Acceptance_v002)  Batch job for runLoopClosureFftnoise_v022, NoiseModel 'endMatched' (TODO v009 N-7).
# Submit FROM the RDS src/ folder that holds runLoopClosureFftnoise_v022.m (Slurm runs the job in the submit directory):
#     cd <RDS project>/src
#     mkdir -p logs
#     sbatch --account=<project code from my_bluebear> runV022HPC_v002.sh "Cook CTRL"
#     sbatch --account=<code> runV022HPC_v002.sh "Hickman PLAC" "Hickman HALO" "Cook CTRL" "Cook ASD"
#     sbatch --account=<code> runV022HPC_v002.sh FraserN48
# Acceptance first (small grid, about 10 min; ask for fewer cores and time):
#     sbatch --account=<code> --cpus-per-task=16 --mem=64G --time=1:00:00 runV022HPC_v002.sh --acceptance "Cook CTRL"
# Watch:   squeue -u $USER        tail -f logs/slurm_v022_<jobid>.out        sacct -j <jobid>        seff <jobid>
# Stop:    scancel <jobid>
# Resume:  a run that PAUSED at the deadline (or was killed) keeps its finished chunks in *_chunks/; submit the same command
#          again and it carries on. Chained: sbatch --dependency=afterany:<jobid> --account=<code> runV022HPC_v002.sh FraserN48
# Settings live in runV022HPC_v002.m (production template, N_BETA 25, N_REPS 200, seed 1729, chunks of 64 trials).
# The account is NOT guessed: pass --account on the sbatch line (or add an #SBATCH --account line once known).
set -euo pipefail

module purge; module load bluebear                    # required first by BlueBEAR
module load bear-apps/2023a
module load MATLAB/2024b-r6

ACCEPT=false
if [[ "${1:-}" == "--acceptance" ]]; then ACCEPT=true; shift; fi
if [[ $# -lt 1 ]]; then echo "usage: sbatch runV022HPC_v002.sh [--acceptance] <dataset> [<dataset> ...]" >&2; exit 2; fi

# per-job temporary folder for the parallel pool's files (BEAR Apps MATLAB page)
export TMPDIR="${TMPDIR:-/tmp}/matlab_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

# MATLAB string array of the datasets, e.g. ["Cook CTRL","Cook ASD"]
DS=""; for d in "$@"; do DS+="\"$d\","; done; DS="[${DS%,}]"

# the job's own end time, so the runner pauses cleanly before Slurm kills it
END=$(squeue -h -j "$SLURM_JOB_ID" -o %e)

echo "Job $SLURM_JOB_ID on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cores, datasets $DS, ends $END, acceptance $ACCEPT"
if $ACCEPT; then
    matlab -batch "runV022HPC_v002($DS, \"\", Acceptance=true)"
else
    matlab -batch "runV022HPC_v002($DS, \"$END\")"
fi
echo "MATLAB exited with status $? at $(date)"
