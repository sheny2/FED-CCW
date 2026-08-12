#!/bin/bash
#SBATCH --job-name=sim_ccw_4_3_10sites_fixed
#SBATCH --array=1-27              # 3 tau x 3 beta_trt x 3 confounding levels
#SBATCH --cpus-per-task=5         # cores per array task (used by mclapply)
#SBATCH --mem=20G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

set -eo pipefail

# --- environment ---------------------------------------------------------
# Load R however your cluster provides it; edit to match your site.
# JHPCE's Conda activation hooks inspect variables that may be unset, so enable
# Bash nounset only after the module has finished loading.
module load R/4.3
set -u

# Keep BLAS/OpenMP from oversubscribing: parallelism is via mclapply only.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Task ${SLURM_ARRAY_TASK_ID} on $(hostname) with ${SLURM_CPUS_PER_TASK} cores"

# --- run one parameter combination ---------------------------------------
Rscript run_sim_array.R "${SLURM_ARRAY_TASK_ID}"

echo "Task ${SLURM_ARRAY_TASK_ID} done."
