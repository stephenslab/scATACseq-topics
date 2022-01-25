#!/bin/bash

#SBATCH --job-name=motif_analysis
#SBATCH --output=motif_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

# This script allocates computing resources (CPUs, memory), loads R,
# and runs postfit_motif_analysis.R.

# Get the command-line arguments.
DAFILE=${1}
GENOME=${2}
SELECTMETHOD=${3}
OUTDIR=${4}

echo dafile=${DAFILE}
echo genome=${GENOME}
echo selectmethod=${SELECTMETHOD}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# path for HOMER
HOMERPATH=/project2/xinhe/software/homer/bin/findMotifsGenome.pl

# Run the R script.
echo "Running postfit_motif_analysis.R"
mkdir -p ${OUTDIR}

Rscript ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.R \
  --DAfile ${DAFILE} \
  --genome ${GENOME} \
  --selectmethod ${SELECTMETHOD} \
  --homerpath ${HOMERPATH} \
  --nc 8 \
  -o ${OUTDIR} > ${OUTDIR}.Rout

