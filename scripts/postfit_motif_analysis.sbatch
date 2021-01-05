#!/bin/bash

#SBATCH --job-name=motif_analysis
#SBATCH --output=motif_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=10:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs diffcount_motif_analysis.R.

# Get the command-line arguments.
DATFILE=${1}
FITFILE=${2}
GENOME=${3}
SELECTMETHOD=${4}
OUTDIR=${5}
echo datfile=${DATFILE}
echo fitfile=${FITFILE}
echo genome=${GENOME}
echo selectmethod=${SELECTMETHOD}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running diffcount_motif_analysis.R"
mkdir -p ${OUTDIR}

Rscript ~/projects/scATACseq-topics/scripts/diffcount_motif_analysis.R \
  --counts ${DATFILE} --modelfit ${FITFILE} \
  --genome ${GENOME} --selectmethod ${SELECTMETHOD} \
  --homerpath /project2/xinhe/software/homer/bin/findMotifsGenome.pl \
  --nc 8 -o ${OUTDIR} > ${OUTDIR}.Rout

