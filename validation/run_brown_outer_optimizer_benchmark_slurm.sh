#!/bin/bash
#SBATCH --partition=sphinx
#SBATCH --account=nlp
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=cdrgam-brown-outer-bench
#SBATCH --output=validation/output/brown_outer_optimizer_benchmark_%j.log

set -euo pipefail

repository="/juice6/u/nlp/climblab/code/cdrgam"
cd "$repository"
export R_ENVIRON_USER=/dev/null
export R_LIBS_USER="$repository/.prototype-lib:$repository/.r-lib"
threads="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS="$threads"
export MKL_NUM_THREADS="$threads"

exec /usr/bin/time -v Rscript validation/brown_outer_optimizer_benchmark.R
