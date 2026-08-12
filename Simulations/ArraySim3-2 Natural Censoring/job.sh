#!/bin/bash
#SBATCH --job-name=sim_ccw
#SBATCH --array=1-27
#SBATCH --cpus-per-task=12         # cores per array task (used by mclapply)
#SBATCH --mem=50G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

# set -eo pipefail

# mkdir logs results

# --- environment ---------------------------------------------------------
# Load R however your cluster provides it; edit to match your site.
module load R/4.3

# Keep BLAS/OpenMP from oversubscribing: parallelism is via mclapply only.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Task ${SLURM_ARRAY_TASK_ID} on $(hostname) with ${SLURM_CPUS_PER_TASK} cores"

# --- run one parameter combination ---------------------------------------
Rscript run_sim_array.R "${SLURM_ARRAY_TASK_ID}"

echo "Task ${SLURM_ARRAY_TASK_ID} done."
