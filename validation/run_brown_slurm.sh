#!/bin/bash
#SBATCH --partition=sphinx
#SBATCH --account=nlp
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=cdrgam-brown
#SBATCH --output=validation/output/brown_%x_%j.log

set -euo pipefail

model_name="${1:?usage: run_brown_slurm.sh MODEL}"
repository="/juice6/u/nlp/climblab/code/cdrgam"
cd "$repository"
export R_LIBS_USER="$repository/.r-lib"
export CDRGAM_BROWN_MODELS="$model_name"
blas_threads="${CDRGAM_BLAS_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="$blas_threads"
export OPENBLAS_NUM_THREADS="$blas_threads"
export MKL_NUM_THREADS="$blas_threads"

exec /usr/bin/time -v Rscript validation/brown_integration.R
