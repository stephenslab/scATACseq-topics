#!/bin/bash

#SBATCH --job-name=cisTopic
#SBATCH --output=cisTopic_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --time=36:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs fit_cisTopic_topics.R.

# Get the command-line arguments.
DATFILE=${1}
NUMITER=${2}
OUTFILE=${3}
echo datfile=${DATFILE}
echo numiter=${NUMITER}
echo outfile=${OUTFILE}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

# Run the R script.
echo "Running fit_cisTopic_topics.R"
export OPENBLAS_NUM_THREADS=1
Rscript ~/projects/scATACseq-topics/scripts/fit_cisTopic_topics.R \
  --counts ${DATFILE} -n ${NUMITER} \
  --nc 10 -o ${OUTFILE}.rds > ${OUTFILE}.Rout
