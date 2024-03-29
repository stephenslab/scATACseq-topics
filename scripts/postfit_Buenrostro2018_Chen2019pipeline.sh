#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# postfit analysis for Buenrostro 2018 single-cell ATAC-seq data sets

# --------------------------------------
# Data processed using Chen2019 pipeline
# --------------------------------------
# Settings
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline/
FIT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit/

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit


# MOTIF ANALYSIS
# --------------
## Compute motif enrichment for each topic using HOMER. Select regions by quantile > 0.99
POSTFIT_DIR=${OUT_DIR}/motifanalysis-Buenrostro2018-k=11-quantile

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       hg19 quantile ${POSTFIT_DIR}

# GENE ANALYSIS
# -------------
## compute gene scores using TSS based method, use original region z-scores, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-none-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS none l2 ${POSTFIT_DIR}

## compute gene scores using TSS based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-abs-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS abs l2 ${POSTFIT_DIR}

## compute gene scores using genebody based method, use original region z-scores, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-abs-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody none l2 ${POSTFIT_DIR}

## compute gene scores using genebody based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-abs-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody abs l2 ${POSTFIT_DIR}


