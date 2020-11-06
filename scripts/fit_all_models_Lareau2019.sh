#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for Lareau 2019 single-cell ATAC-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch


# all Lareau_2019_bonemarrow data
# ============================================================================================================================
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/bone_marrow/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow
mkdir -p ${OUT_DIR}

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019/Lareau_2019_bonemarrow
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019/Lareau_2019_bonemarrow

# "Pre-fit" factorizations to the Lareau_2019_bonemarrow data.
# Computation took 33196 seconds (9 hrs) for 300 iterations with k = 3 (using 10 cpus).
# Computation took 43196 seconds (12 hrs) for 300 iterations with k = 8 (using 10 cpus).
# Computation took 55288 seconds (15 hrs) for 300 iterations with k = 13 (using 10 cpus).
#                                 data                                     k  numiter outfile
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  2  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  3  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  4  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  5  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  6  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  7  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  8  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  9  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  10 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  11 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  12 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  13 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13

# Fit factorizations to Lareau_2019_bonemarrow data, with and without extrapolation.
# Computation took 124272 seconds (34 hrs) for 600 iterations with k = 13, and scd method (using 10 cpus).
#                    data                                    prefitfile                                  k method numiter ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=2

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=3 3 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=3

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=4 4 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=4

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=5 5 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=5

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=6 6 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=6

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=7 7 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=7

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=8 8 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=8

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=9 9 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=9

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=10 10 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=10

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=11 11 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=11

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=12 12 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=12

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=13 13 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Lareau2019_bonemarrow

# Lareau_2019_bonemarrow_resting data
# ============================================================================================================================

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/bone_marrow/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow_resting
mkdir -p ${OUT_DIR}

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019_bonemarrow_resting
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019_bonemarrow_resting

# "Pre-fit" factorizations to the Lareau_2019_bonemarrow_resting data.
# Computation took 21000 seconds (6 hrs) for 300 iterations with k = 13 (using 10 cpus).
#                                 data                                             k  numiter outfile
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  2  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  3  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  4  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  5  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  6  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  7  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  8  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  9  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  10 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  11 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  12 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData  13 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13

# Fit factorizations to Lareau_2019_bonemarrow_resting data, with and without extrapolation.
#                    data                                            prefitfile                                          k method numiter ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=2 2 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=2

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=3 3 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=3

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=4 4 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=4

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=5 5 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=5

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=6 6 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=6

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=7 7 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=7

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=8 8 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=8

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=9 9 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=9

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=10 10 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=10

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=11 11 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=11

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=12 12 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=12

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-em-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-ccd-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_resting.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_resting-k=13 13 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_resting-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow_resting
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Lareau2019_bonemarrow_resting

# Lareau_2019_bonemarrow_stimulated data
# ============================================================================================================================

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/bone_marrow/processed_data
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow_stimulated
mkdir -p ${OUT_DIR}

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019_bonemarrow_stimulated
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Lareau_2019_bonemarrow_stimulated

# "Pre-fit" factorizations to the Lareau_2019_bonemarrow_stimulated data.
# Computation took 34000 seconds (9 hrs) for 300 iterations with k = 13 (using 10 cpus).
#                                 data                                                k  numiter outfile
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  2  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  3  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  4  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  5  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  6  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  7  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  8  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  9  300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  10 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  11 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  12 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData  13 300     ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13

# Fit factorizations to Lareau_2019_bonemarrow_stimulated data, with and without extrapolation.
#                    data                                               prefitfile                                             k method numiter ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=2 2 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=2

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=3
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=3 3 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=3

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=4
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=4 4 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=4

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=5
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=5 5 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=5

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=6
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=6 6 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=6

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=7
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=7 7 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=7

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=8
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=8 8 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=8

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=9
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=9 9 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=9

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=10
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=10 10 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=10

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=11
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=11 11 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=11

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=12
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=12 12 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=12

sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 em     600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 ccd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 scd    600     no ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 em     600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-em-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 ccd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-ccd-ex-k=13
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow_stimulated.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow_stimulated-k=13 13 scd    600     yes ${OUT_DIR}/fit-Lareau2019_bonemarrow_stimulated-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019/Lareau_2019_bonemarrow_stimulated
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Lareau2019_bonemarrow_stimulated

# # Lareau_2019_mousebrain data
# DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/mouse_brain/processed_data
# OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Lareau_2019
# NUMITER=10
# mkdir -p ${OUT_DIR}
#
# cd /project2/mstephens/kevinluo/scATACseq-topics/log
#
# # "Pre-fit" factorizations to the Lareau_2019_mousebrain data.
# #                       data                                     k numiter    outfile
# sbatch ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData  2 ${NUMITER} ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2
#
# # Fit factorizations to Lareau_2019_mousebrain data, with and without extrapolation.
# #                    data                                    prefitfile                                  k method numiter     ex outfile
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 em     ${NUMITER}  no fit-Lareau2019_mousebrain-em-k=2
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 ccd    ${NUMITER}  no fit-Lareau2019_mousebrain-ccd-k=2
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 scd    ${NUMITER}  no fit-Lareau2019_mousebrain-scd-k=2
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 em     ${NUMITER} yes fit-Lareau2019_mousebrain-em-ex-k=2
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 ccd    ${NUMITER} yes fit-Lareau2019_mousebrain-ccd-ex-k=2
# sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 scd    ${NUMITER} yes fit-Lareau2019_mousebrain-scd-ex-k=2

