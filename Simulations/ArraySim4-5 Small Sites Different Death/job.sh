#!/bin/bash
#SBATCH --job-name=sim_ccw_4_5
#SBATCH --array=1-12
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

set -eo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

module load R/4.3
set -u

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Task ${SLURM_ARRAY_TASK_ID} on $(hostname) with ${SLURM_CPUS_PER_TASK} cores"
Rscript run_sim_array.R "${SLURM_ARRAY_TASK_ID}"
echo "Task ${SLURM_ARRAY_TASK_ID} done."
