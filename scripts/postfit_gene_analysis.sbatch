#!/bin/bash

#SBATCH --job-name=gene_analysis
#SBATCH --output=gene_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

# Get the command-line arguments.
DATFILE=${1}
FITFILE=${2}
GENOME=${3}
GENESCOREMETHOD=${4}
TRANSFORM=${5}
NORMALIZATION=${6}
OUTDIR=${7}

echo datfile=${DATFILE}
echo fitfile=${FITFILE}
echo genome=${GENOME}
echo genescoremethod=${GENESCOREMETHOD}
echo transform=${TRANSFORM}
echo normalization=${NORMALIZATION}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running diffcount_genescore_gsea.R"
mkdir -p ${OUTDIR}

Rscript ~/projects/scATACseq-topics/scripts/diffcount_genescore.R \
  --counts ${DATFILE} \
  --modelfit ${FITFILE} \
  --genome ${GENOME} \
  --genescoremethod ${GENESCOREMETHOD} \
  --transform ${TRANSFORM} \
  --normalization ${NORMALIZATION} \
  --nc 8 -o ${OUTDIR} > ${OUTDIR}.genescore.Rout

Rscript ~/projects/scATACseq-topics/scripts/perform_gsea.R \
  --genescore ${OUTDIR}/genescore_result.rds \
  --genome ${GENOME} \
  --nc 8 -o ${OUTDIR} > ${OUTDIR}.gsea.Rout
