#!/bin/bash

#SBATCH --job-name=fitNMF
#SBATCH --output=fitNMF_%J.out
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G

# This script allocates computing resources (CPUs, memory), loads R,
# and runs prefit and fit poisson NMF.

# Get the command-line arguments.
K=${1}
NUMPREFIT=300
NUMITER=600
echo k=${K}
echo numprefit=${NUMPREFIT}
echo numiter=${NUMITER}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# Run the R script.
echo "Running fit_poisson_nmf_Buenrostro2018_filtered.R"
export OPENBLAS_NUM_THREADS=1
Rscript ~/projects/scATACseq-topics/scripts/fit_poisson_nmf_Buenrostro2018_filtered.R \
  -k ${K} --numprefit ${NUMPREFIT} -n ${NUMITER} \
  --nc 10 > ${OUTFILE}.Rout
