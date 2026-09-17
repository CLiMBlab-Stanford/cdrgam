#!/bin/bash
#SBATCH --partition=sphinx
#SBATCH --account=nlp
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=cdrgam-brown-trust
#SBATCH --output=validation/output/brown_%x_%j.log

set -euo pipefail

model_name="${1:?usage: run_brown_trust_slurm.sh MODEL}"
repository="/juice6/u/nlp/climblab/code/cdrgam"
cd "$repository"

export R_ENVIRON_USER=/dev/null
trust_library="${CDRGAM_TRUST_LIBRARY:-$repository/.trust-lib}"
export R_LIBS_USER="$trust_library:$repository/.r-lib"
export CDRGAM_BROWN_MODELS="$model_name"
export CDRGAM_BROWN_FULL=1
export CDRGAM_BROWN_BACKEND=sparse_trust
export CDRGAM_BROWN_RUN_LABEL="${CDRGAM_BROWN_RUN_LABEL:-full_subject_irf_trust}"
export CDRGAM_BROWN_SAVE_FITS=1
export CDRGAM_BROWN_CHECKPOINT=1
export CDRGAM_BROWN_SOLVER_TRACE=1
export CDRGAM_BROWN_OPTIMIZER_MAXIT=300
export CDRGAM_BROWN_OPTIMIZER_GRADIENT_TOLERANCE=2e-4

blas_threads="${CDRGAM_BLAS_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="$blas_threads"
export OPENBLAS_NUM_THREADS="$blas_threads"
export MKL_NUM_THREADS="$blas_threads"

exec /usr/bin/time -v Rscript validation/brown_integration.R
