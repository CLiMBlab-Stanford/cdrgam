#!/bin/bash
#SBATCH --partition=sphinx
#SBATCH --account=nlp
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=cdrgam-hessian-bench
#SBATCH --output=validation/output/hessian_benchmark_%j.log

set -euo pipefail

repository="/juice6/u/nlp/climblab/code/cdrgam"
cd "$repository"
export R_LIBS_USER="$repository/.r-lib"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

exec /usr/bin/time -v Rscript validation/hessian_benchmark.R
