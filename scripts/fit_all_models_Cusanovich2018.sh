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
# Computation took 64072 seconds (17 hrs) for 200 iterations with k = 3 (using 10 cpus).
# Computation took 79812 seconds (22 hrs) for 200 iterations with k = 8 (using 10 cpus).
# Computation took 103200 seconds (28 hrs) for 200 iterations with k = 13 (using 10 cpus).
#                                 data                              k  numiter outfile
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=45G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
#                              data                             prefitfile                           k method numiter ex  outfile
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
