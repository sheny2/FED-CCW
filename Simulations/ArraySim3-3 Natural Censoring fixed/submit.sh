#!/bin/bash
set -euo pipefail

# SLURM opens the log files before job.sh starts, so these directories must
# exist before sbatch is called.
mkdir -p logs results
sbatch job.sh
