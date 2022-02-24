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
# Computation took 148200 seconds (41 hrs) for 200 iterations with k = 25 (using 10 cpus).
# Computation took 150592 seconds (42 hrs) for 200 iterations with k = 30 (using 10 cpus).
# Computation took 207743 seconds (58 hrs) for 200 iterations with k = 40 (using 28 cpus).

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
sbatch --mem=55G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  14 200     ${OUT_DIR}/prefit-Cusanovich2018-k=14
sbatch --mem=50G ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  15 200     ${OUT_DIR}/prefit-Cusanovich2018-k=15

sbatch --mem=70G --partition=gilad --account=pi-gilad ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  20 200     ${OUT_DIR}/prefit-Cusanovich2018-k=20
sbatch --mem=70G --partition=xinhe --account=pi-xinhe ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  25 200     ${OUT_DIR}/prefit-Cusanovich2018-k=25
sbatch --mem=70G --partition=mstephens ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  30 200     ${OUT_DIR}/prefit-Cusanovich2018-k=30
sbatch --mem=100G --partition=mstephens --cpus-per-task=28 --exclusive ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  35 200     ${OUT_DIR}/prefit-Cusanovich2018-k=35
sbatch --mem=100G --partition=gilad --account=pi-gilad --cpus-per-task=28 ${SCRIPT_PREFIT} ${DAT_DIR}/Cusanovich_2018.RData  40 200     ${OUT_DIR}/prefit-Cusanovich2018-k=40

# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for 250 iterations with k = 7
# Computation took 228004 seconds (63 hrs) for 300 iterations with k = 25
# Computation took 229813 seconds (64 hrs) for 300 iterations with k = 30 (using 10 cpus)
# Computation took 282852 seconds (78 hrs) for 300 iterations with k = 35 (using 28 cpus)
# Computation took 322392 seconds (89 hrs) for 300 iterations with k = 40 (using 28 cpus)

#                              data                             prefitfile                           k method numiter ex  outfile
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=2
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=2
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=2
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=2
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=2
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=3
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=3
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=3
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=3
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=3
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=4
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=4
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=4
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=4
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=4
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=5
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=6
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=6
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=6
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=6
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=6
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=7
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=8
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=8
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=8
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=8
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=8
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=9
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=9
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=9
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=9
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=9
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9

sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=10
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=10
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=10
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=10
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=10
sbatch --mem=50G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10

sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=11
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=11
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=11
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=11
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=11
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11

sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=12
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=12
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=12
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=12
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=12
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12

sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=13
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=13
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=13
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=13
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=13
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=14
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=14
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=14
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=14
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=14
sbatch --mem=55G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=14 14 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=14
#
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 em     250     no  ${OUT_DIR}/fit-Cusanovich2018-em-k=15
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 ccd    250     no  ${OUT_DIR}/fit-Cusanovich2018-ccd-k=15
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 scd    250     no  ${OUT_DIR}/fit-Cusanovich2018-scd-k=15
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 em     250     yes ${OUT_DIR}/fit-Cusanovich2018-em-ex-k=15
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 ccd    250     yes ${OUT_DIR}/fit-Cusanovich2018-ccd-ex-k=15
sbatch --mem=56G ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=15 15 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=15

sbatch --mem=80G --partition=mstephens ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=20 20 scd    300     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=20
sbatch --mem=90G --partition=gilad --account=pi-gilad ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=25 25 scd    300     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=25
sbatch --mem=90G --partition=mstephens ${SCRIPT_FIT} ${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=30 30 scd    300     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=30

sbatch --mem=100G --partition=mstephens --cpus-per-task=28 --exclusive ${SCRIPT_FIT} \
${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=35 35 scd    300     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=35

sbatch --mem=100G --partition=mstephens --cpus-per-task=28 --exclusive ${SCRIPT_FIT} \
${DAT_DIR}/Cusanovich_2018.RData ${OUT_DIR}/prefit-Cusanovich2018-k=40 40 scd    300     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=40

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018
