#!/bin/bash

#SBATCH --job-name=DA_analysis
#SBATCH --output=DA_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G

# This script allocates computing resources (CPUs, memory), loads R,
# and runs DA_analysis.R.

# Get the command-line arguments.
DATFILE=${1}
FITFILE=${2}
NS=${3}
NSPLIT=${4}
OUTDIR=${5}
echo datfile=${DATFILE}
echo fitfile=${FITFILE}
echo ns=${NS}
echo nsplit=${NSPLIT}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# Run the R script.
echo "Running DA_analysis.R"
mkdir -p ${OUTDIR}

Rscript ~/projects/scATACseq-topics/scripts/DA_analysis.R \
  --counts ${DATFILE} --modelfit ${FITFILE} \
  --nc 20 \
  --ns ${NS} \
  --nsplit ${NSPLIT} \
  -o ${OUTDIR} > ${OUTDIR}.DA.${NS}iters.Rout

