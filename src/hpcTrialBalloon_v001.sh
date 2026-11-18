#!/bin/bash
#SBATCH --job-name=balloon_v022
#SBATCH --qos=bbshort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --output=logs/slurm_balloon_%j.out
#
# hpcTrialBalloon_v001.sh  (2026-11-15)  Ten-minute test of the BlueBEAR batch route (docs/README_BlueBEARBatch_v001.md).
# Same module lines, $TMPDIR handling and deadline logic as runV022HPC_v001.sh; runs hpcTrialBalloon_v001.m, which checks
# MATLAB in batch, the licence, the BEAR-style pool, every file the v022 run needs, the Cook CTRL import, the end-matched
# pieces on one real trial against the M5's values, deadline parsing and write access.
# Submit from the RDS src/ folder:
#     mkdir -p logs
#     sbatch --account=<project code from my_bluebear> hpcTrialBalloon_v001.sh
# Then:  squeue -u $USER      cat logs/slurm_balloon_<jobid>.out      sacct -j <jobid>
# Pass:  the .out file ends with "TRIAL BALLOON PASSED" and sacct shows COMPLETED.
set -euo pipefail

module purge; module load bluebear
module load bear-apps/2023a
module load MATLAB/2024b-r6

export TMPDIR="${TMPDIR:-/tmp}/matlab_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
END=$(squeue -h -j "$SLURM_JOB_ID" -o %e)

echo "Balloon job $SLURM_JOB_ID on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cores, ends $END, TMPDIR $TMPDIR"
matlab -batch "hpcTrialBalloon_v001(\"$END\")"
echo "MATLAB exited with status $? at $(date)"
