#!/bin/bash
#SBATCH --partition=sphinx
#SBATCH --account=nlp
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --job-name=cdrgam-effect-calibration
#SBATCH --output=validation/output/effect_calibration_%j.log
#
# Usage: sbatch validation/run_effect_calibration_slurm.sh [replicates] [first_seed] [backend]
# The current checkout is installed into a job-local library, so the run
# tests the code as checked out and leaves shared libraries alone.

set -euo pipefail

repository="/juice6/u/nlp/climblab/code/cdrgam"
cd "$repository"
environment_bin="/nlp/scr/cashain/miniconda3/envs/cdrgam/bin"
export PATH="$environment_bin:$PATH"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

library=$(mktemp -d "${TMPDIR:-/tmp}/cdrgam-calibration-lib.XXXXXX")
trap 'rm -rf "$library"' EXIT
"$environment_bin/R" CMD INSTALL --no-test-load -l "$library" . >/dev/null
export R_LIBS="$library"

exec /usr/bin/time -v "$environment_bin/Rscript" validation/effect_calibration.R \
    "${1:-500}" "${2:-1}" "${3:-sparse}"
