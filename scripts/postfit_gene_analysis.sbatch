#!/bin/bash

#SBATCH --job-name=gene_analysis
#SBATCH --output=gene_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=10:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs diffcount_genescore_gsea.R.

# Get the command-line arguments.
DATFILE=${1}
FITFILE=${2}
GENESET=${3}
GENOME=${4}
GENESCOREMETHOD=${5}
NORMALIZATION=${6}
OUTDIR=${7}
echo datfile=${DATFILE}
echo fitfile=${FITFILE}
echo geneset=${GENESET}
echo genome=${GENOME}
echo genescoremethod=${GENESCOREMETHOD}
echo normalization=${NORMALIZATION}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running diffcount_genescore_gsea.R"
Rscript ~/projects/scATACseq-topics/scripts/diffcount_genescore_gsea.R \
  --counts ${DATFILE} --modelfit ${FITFILE} \
  --geneset ${GENESET} --genome ${GENOME} \
  --genescoremethod {GENESCOREMETHOD} --normalization {NORMALIZATION} \
  --nc 8 -o ${OUTDIR} > ${OUTDIR}.Rout
