#!/bin/bash

#SBATCH --job-name=fitNMF
#SBATCH --output=fitNMF_%J.out
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G

# This script allocates computing resources (CPUs, memory), loads R,
# and runs prefit and fit poisson NMF.

# Get the command-line arguments.
TISSUE=${1}
K=${2}

DATFILE=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data/Cusanovich_2018_${TISSUE}.RData
NUMPREFIT=300
NUMITER=600

OUTDIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018/tissues
mkdir -p ${OUTDIR}
ROUTFILE=${OUTDIR}/fit_poisson_nmf_Cusanovich2018_${TISSUE}.Rout

echo k=${K}
echo tissue=${TISSUE}
echo numprefit=${NUMPREFIT}
echo numiter=${NUMITER}
echo logfile=${ROUTFILE}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# Run the R script.
echo "Running fit_poisson_nmf_Cusanovich2018_tissue.R"
export OPENBLAS_NUM_THREADS=1
Rscript ~/projects/scATACseq-topics/scripts/fit_poisson_nmf_Cusanovich2018_tissue.R \
  --counts ${DATFILE} --tissue ${TISSUE} -k ${K} \
  --numprefit ${NUMPREFIT} -n ${NUMITER} --nc 10 \
  > ${ROUTFILE}
