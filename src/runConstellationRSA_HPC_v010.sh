#!/bin/bash
#SBATCH --job-name=rsa010
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/slurm_rsa010_%j.out
#SBATCH --mail-type=END,FAIL
#
# runConstellationRSA_HPC_v010.sh (2026-11-16)  Cross-trial RSA (Mantel) on the v028 end-matched corpus, full N, NPerm 10000.
# Submit FROM the run folder rsa_v010/ (RDS _N48 tree; it holds this script, runConstellationRSA_HPC_v010.m, buildConstellationRDM_v002.m,
# functions/loadV028_v001.m, mantelToolbox/ and symlinks to the loopClosureResults_*_v028.mat files), always with --account=fraserds-mpo-evaluation:
#   sbatch --account=fraserds-mpo-evaluation runConstellationRSA_HPC_v010.sh                              primary arm, window basis
#   sbatch --account=fraserds-mpo-evaluation runConstellationRSA_HPC_v010.sh --obs whole                  whole-trial beta_obs (sensitivity)
#   sbatch --account=fraserds-mpo-evaluation runConstellationRSA_HPC_v010.sh --arm shapeXu                old-noise arm (sensitivity)
# Outputs constellationRSA_HPC_v010_<tag>.mat beside the script. Watch: squeue -u $USER; tail -f logs/slurm_rsa010_<jobid>.out
set -euo pipefail
module purge; module load bluebear
module load bear-apps/2023a
module load MATLAB/2024b-r6
BASIS=v028; ARM=primary; OBS=window
while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
        --basis) BASIS="$2"; shift 2 ;;
        --arm) ARM="$2"; shift 2 ;;
        --obs) OBS="$2"; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
done
mkdir -p logs
export TMPDIR="${TMPDIR:-/tmp}/matlab_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
echo "Job $SLURM_JOB_ID on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cores, basis $BASIS arm $ARM obs $OBS, start $(date)"
matlab -batch "RSA_BASIS='$BASIS'; RSA_ARM='$ARM'; RSA_OBS='$OBS'; runConstellationRSA_HPC_v010"
echo "MATLAB finished at $(date)"
