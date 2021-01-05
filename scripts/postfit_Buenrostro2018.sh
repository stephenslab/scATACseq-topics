#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# postfit analysis for Buenrostro 2018 single-cell ATAC-seq data sets

# GENE ANALYSIS
# -------------

# Data processed using Chen2019 pipeline
# --------------------------------------
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit

## compute gene scores using genebody based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-l2

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody l2 ${POSTFIT_DIR}

## compute gene scores using genebody based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-sum

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody sum ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-l2

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS l2 ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-sum

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS sum ${POSTFIT_DIR}

# Data processed using chromVAR scPeaks
# --------------------------------------
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline/chromVAR/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_chromVAR_scPeaks/binarized/

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit

## compute gene scores using genebody based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-l2

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody l2 ${POSTFIT_DIR}

## compute gene scores using genebody based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-genebody-sum

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 genebody sum ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the l2 norm of the weights
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-l2

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS l2 ${POSTFIT_DIR}

## compute gene scores using TSS based method, and normalize by the sum of weights.
POSTFIT_DIR=${OUT_DIR}/geneanalysis-Buenrostro2018-k=11-TSS-sum

sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_gene_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       /project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_human.RData \
       hg19 TSS sum ${POSTFIT_DIR}

# MOTIF ANALYSIS
# --------------
# Data processed using Chen2019 pipeline
# --------------------------------------
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit

## Compute motif enrichment for each topic using HOMER. Select regions by quantile > 0.99
sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       hg19 quantile ${OUT_DIR}/motifanalysis-Buenrostro2018-k=11-quantile

# Data processed using chromVAR scPeaks
# --------------------------------------
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline/chromVAR/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_chromVAR_scPeaks/binarized/

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018/postfit

## Compute motif enrichment for each topic using HOMER. Select regions by quantile > 0.99
sbatch --mem=20G ~/projects/scATACseq-topics/scripts/postfit_motif_analysis.sbatch \
       ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData \
       ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11.rds \
       hg19 quantile ${OUT_DIR}/motifanalysis-Buenrostro2018-k=11-quantile