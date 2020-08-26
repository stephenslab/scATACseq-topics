#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for Cusanovich 2018 single-cell ATAC-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018
mkdir -p ${OUT_DIR}

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

# "Pre-fit" factorizations to the Cusanovich_2018 data.
#                                 data                              k  numiter outfile
sbatch --mem=40G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  2  100     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  3  100     ${OUT_DIR}/prefit-Cusanovich2018-k=3

sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5

sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  6  300     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  7  300     ${OUT_DIR}/prefit-Cusanovich2018-k=7


# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
#                    data                             prefitfile                           k method numiter     ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     500         no ${OUT_DIR}/fit-Cusanovich2018-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    500         no ${OUT_DIR}/fit-Cusanovich2018-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    500         no ${OUT_DIR}/fit-Cusanovich2018-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     500         yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    500         yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    500         yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2


