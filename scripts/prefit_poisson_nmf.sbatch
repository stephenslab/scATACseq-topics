#!/bin/bash

#SBATCH --job-name=prefitNMF
#SBATCH --output=prefitNMF_%J.out
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G

# This script allocates computing resources (CPUs, memory), loads R,
# and runs prefit_poisson_nmf.R. See fit_all_models.sh for examples
# illustrating how to use this script.

# Get the command-line arguments.
DATFILE=${1}
K=${2}
NUMITER=${3}
OUTFILE=${4}
echo datfile=${DATFILE}
echo k=${K}
echo numiter=${NUMITER}
echo outfile=${OUTFILE}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running prefit_poisson_nmf.R"
export OPENBLAS_NUM_THREADS=1
Rscript ~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.R \
  --counts ${DATFILE} -k ${K} -n ${NUMITER} \
  --nc 10 -o ${OUTFILE}.rds > ${OUTFILE}.Rout
