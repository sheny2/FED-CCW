#!/bin/bash
set -euo pipefail

# SLURM opens the log files before job.sh starts, so these directories must
# exist before sbatch is called.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p logs results
sbatch job.sh
