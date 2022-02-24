#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# postfit analysis for Cusanovich 2018 single-cell ATAC-seq data sets

# Settings
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
FIT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018/postfit_v2

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit_v2
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/postfit_v2

# DIFF ACCESSIBILITY (DA) ANALYSIS
# ----------------------------------------
## Compute differential accessibility across regions
DA_DIR=${OUT_DIR}/DAanalysis-Cusanovich2018-k=13

sbatch --mem=56G ~/projects/scATACseq-topics/scripts/postfit_DA_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${FIT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       10000 100 ${DA_DIR}

sbatch --mem=56G ~/projects/scATACseq-topics/scripts/postfit_DA_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${FIT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       1000 100 ${DA_DIR}

sbatch --mem=56G ~/projects/scATACseq-topics/scripts/postfit_DA_analysis.sbatch \
       ${DAT_DIR}/Cusanovich_2018.RData \
       ${FIT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
       5000 100 ${DA_DIR}

# sinteractive -p broadwl -c 20 --account=pi-mstephens --mem=40G --time=18:0:0
Rscript ~/projects/scATACseq-topics/scripts/DA_analysis.R \
  --counts ${DAT_DIR}/Cusanovich_2018.RData \
  --modelfit ${FIT_DIR}/fit-Cusanovich2018-scd-ex-k=13.rds \
  --nc 20 \
  --ns 100 \
  --nsplit 100 \
  -o ${DA_DIR}

# MOTIF ANALYSIS
# --------------

## Compute motif enrichment for each topic using HOMER. Select regions by quantile
MOTIFANALYSIS_DIR=${OUT_DIR}/motifanalysis-Cusanovich2018-k=13-quantile

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 quantile ${MOTIFANALYSIS_DIR}

## Compute motif enrichment for each topic using HOMER. Select regions by zscore
MOTIFANALYSIS_DIR=${OUT_DIR}/motifanalysis-Cusanovich2018-k=13-zscore

sbatch --mem=40G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 zscore ${MOTIFANALYSIS_DIR}


# GENE ANALYSIS
# --------------

## compute gene scores using TSS based method, use original region z-scores, and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-Z-l2
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 TSS none l2 ${GENEANALYSIS_DIR}

## compute gene scores using TSS based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-absZ-l2
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 TSS abs l2 ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, use original region z-scores, and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-Z-l2
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 genebody none l2 ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-absZ-l2
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 genebody abs l2 ${GENEANALYSIS_DIR}



## compute gene scores using TSS based method, use original region z-scores, and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-Z-sum
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 TSS none sum ${GENEANALYSIS_DIR}

## compute gene scores using TSS based method, transform region z-scores using abs(z), and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-TSS-absZ-sum
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 TSS abs sum ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, use original region z-scores, and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-Z-sum
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 genebody none sum ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, transform region z-scores using abs(z), and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Cusanovich2018-k=13-genebody-absZ-sum
sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       mm9 genebody abs sum ${GENEANALYSIS_DIR}
