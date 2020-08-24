#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for Cusanovich 2018 single-cell ATAC-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/output/
NUMITER=10

cd /project2/mstephens/kevinluo/scATACseq-topics/data/log
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the Cusanovich_2018 data.
#                       data                              k numiter    outfile
sbatch ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  2 ${NUMITER} ${OUT_DIR}/prefit-Cusanovich2018-k=2

# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
#                    data                             prefitfile                           k method numiter     ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     ${NUMITER}  no ${OUT_DIR}/fit-Cusanovich2018-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    ${NUMITER}  no ${OUT_DIR}/fit-Cusanovich2018-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    ${NUMITER}  no ${OUT_DIR}/fit-Cusanovich2018-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     ${NUMITER} yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    ${NUMITER} yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    ${NUMITER} yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2


