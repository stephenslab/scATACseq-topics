#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# postfit analysis for Buenrostro 2018 single-cell ATAC-seq data sets

# -----------------------------------------------------
# Use data processed from the Buenrostro 2018 paper
# -----------------------------------------------------
# Settings
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data
FIT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit_v2
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit_v2

# K = 10
# -------

# DIFF ACCESSIBILITY (DA) ANALYSIS
# ----------------------------------------
## Compute differential accessibility across regions
DA_DIR=${OUT_DIR}/DAanalysis-Buenrostro2018-k=10

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_DA_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${FIT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=10.rds \
       10000 100 ${DA_DIR}


# MOTIF ANALYSIS
# --------------
## Compute motif enrichment for each topic using HOMER. Select regions by quantile > 0.99
MOTIFANALYSIS_DIR=${OUT_DIR}/motifanalysis-Buenrostro2018-k=10-quantile

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 quantile ${MOTIFANALYSIS_DIR}

## Compute motif enrichment for each topic using HOMER. Select regions by zscore
MOTIFANALYSIS_DIR=${OUT_DIR}/motifanalysis-Buenrostro2018-k=10-zscore

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 zscore ${MOTIFANALYSIS_DIR}

# GENE ANALYSIS
# -------------
## compute gene scores using TSS based method, use original region z-scores, and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-TSS-Z-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 TSS none l2 ${GENEANALYSIS_DIR}

## compute gene scores using TSS based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-TSS-absZ-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 TSS abs l2 ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, use original region z-scores, and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-genebody-Z-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 genebody none l2 ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, transform region z-scores using abs(z), and normalize by the l2 norm of the weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-genebody-absZ-l2

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 genebody abs l2 ${GENEANALYSIS_DIR}


## compute gene scores using TSS based method, use original region z-scores, and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-TSS-Z-sum

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 TSS none sum ${GENEANALYSIS_DIR}

## compute gene scores using TSS based method, transform region z-scores using abs(z), and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-TSS-absZ-sum

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 TSS abs sum ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, use original region z-scores, and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-genebody-Z-sum

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 genebody none sum ${GENEANALYSIS_DIR}

## compute gene scores using genebody based method, transform region z-scores using abs(z), and normalize by the sum of weights
GENEANALYSIS_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=10-genebody-absZ-sum

sbatch --mem=30G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis_v2.sbatch \
       ${DA_DIR}/DA_regions_topics_10000iters.rds \
       hg19 genebody abs sum ${GENEANALYSIS_DIR}
