#!/bin/bash

#SBATCH --job-name=gene_analysis
#SBATCH --output=gene_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

# Get the command-line arguments.
DAFILE=${1}
GENOME=${2}
GENESCOREMETHOD=${3}
TRANSFORM=${4}
NORMALIZATION=${5}
OUTDIR=${6}

echo dafile=${DAFILE}
echo genome=${GENOME}
echo genescoremethod=${GENESCOREMETHOD}
echo transform=${TRANSFORM}
echo normalization=${NORMALIZATION}
echo outdir=${OUTDIR}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# Run the R script.
mkdir -p ${OUTDIR}

echo "Running postfit_genescore.R"

Rscript ~/projects/scATACseq-topics/scripts/postfit_genescore.R \
  --DAfile ${DAFILE} \
  --genome ${GENOME} \
  --genescoremethod ${GENESCOREMETHOD} \
  --transform ${TRANSFORM} \
  --normalization ${NORMALIZATION} \
  -o ${OUTDIR} > ${OUTDIR}.genescore.Rout

echo "Running perform_gsea.R"

Rscript ~/projects/scATACseq-topics/scripts/perform_gsea.R \
  --genescore ${OUTDIR}/genescore_result.rds \
  --genome ${GENOME} \
  --nc 8 \
  -o ${OUTDIR} > ${OUTDIR}.gsea.Rout
