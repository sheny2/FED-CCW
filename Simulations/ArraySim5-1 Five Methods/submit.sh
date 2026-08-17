#!/bin/bash
# Submit one SLURM array task per row of the grid defined in params.R.
set -euo pipefail

n_tasks=$(Rscript -e 'source("params.R"); validate_params(); cat(nrow(simulation_grid()))')
sbatch --array="1-${n_tasks}" job.sh
