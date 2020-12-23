#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# postfit analysis for Cusanovich 2018 single-cell ATAC-seq data sets
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit

## compute gene scores using genebody based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-l2
mkdir -p ${POSTFIT_DIR}
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_mouse.RData \
       mm9 genebody l2 ${POSTFIT_DIR}

## compute gene scores using genebody based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-sum
mkdir -p ${POSTFIT_DIR}
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_mouse.RData \
       mm9 genebody sum ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-l2
mkdir -p ${POSTFIT_DIR}
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_mouse.RData \
       mm9 TSS l2 ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-sum
mkdir -p ${POSTFIT_DIR}
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_mouse.RData \
       mm9 TSS sum ${POSTFIT_DIR}


# MOTIF ANALYSIS
# --------------

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit

## Compute motif enrichment for each topic using HOMER. Select regions by quantile
POSTFIT_DIR=${OUT_DIR}/motifanalysis-Cusanovich2018-k=13-quantile
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       mm9 quantile ${POSTFIT_DIR}

## Compute motif enrichment for each topic using HOMER. Select regions by zscore
POSTFIT_DIR=${OUT_DIR}/motifanalysis-Cusanovich2018-k=13-zscore
mkdir -p ${POSTFIT_DIR}
ln -sf ${OUT_DIR}/diffcount-Cusanovich2018-13topics.rds ${POSTFIT_DIR}/diffcount_regions_topics.rds

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       mm9 zscore ${POSTFIT_DIR}

