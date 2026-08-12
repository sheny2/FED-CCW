#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ ! -f external_reference_numerator.rds ]]; then
  echo "Missing external_reference_numerator.rds" >&2
  echo "Run: Rscript prepare_external_reference.R" >&2
  exit 1
fi

mkdir -p logs results
sbatch job.sh
