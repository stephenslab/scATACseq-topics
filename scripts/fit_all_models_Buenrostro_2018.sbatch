#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for Buenrostro 2018 single-cell ATAC-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data/

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Buenrostro_2018

# Binarized counts
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the Buenrostro_2018 binarized counts data.
#                                 data                                               k  numiter outfile
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  2  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=2
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  3  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=3
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  4  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=4
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  5  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=5
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  6  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=6
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  7  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=7
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  8  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=8
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  9  300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=9
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  10 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=10
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  11 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=11
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  12 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=12
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  13 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=13
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  14 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=14
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData  15 300     ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=15

# Fit factorizations to the Buenrostro_2018 binarized counts data
#                              data                                              prefitfile                                      k  method numiter ex  outfile
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=2  2  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=2
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=3  3  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=3
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=4  4  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=4
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=5  5  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=5
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=6  6  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=6
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=7  7  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=7
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=8  8  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=8
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=9  9  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=9
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=10 10 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=10
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=11 11 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=11
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=12 12 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=12
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=13 13 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=13
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=14 14 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=14
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_binarized_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-binarized-k=15 15 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-binarized-scd-ex-k=15

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Buenrostro2018-binarized

# Unbinarized counts
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/unbinarized/
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the Buenrostro_2018 counts data.
#                                 data                                     k  numiter outfile
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  2  300     ${OUT_DIR}/prefit-Buenrostro2018-k=2
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  3  300     ${OUT_DIR}/prefit-Buenrostro2018-k=3
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  4  300     ${OUT_DIR}/prefit-Buenrostro2018-k=4
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  5  300     ${OUT_DIR}/prefit-Buenrostro2018-k=5
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  6  300     ${OUT_DIR}/prefit-Buenrostro2018-k=6
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  7  300     ${OUT_DIR}/prefit-Buenrostro2018-k=7
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  8  300     ${OUT_DIR}/prefit-Buenrostro2018-k=8
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  9  300     ${OUT_DIR}/prefit-Buenrostro2018-k=9
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  10 300     ${OUT_DIR}/prefit-Buenrostro2018-k=10
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  11 300     ${OUT_DIR}/prefit-Buenrostro2018-k=11
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  12 300     ${OUT_DIR}/prefit-Buenrostro2018-k=12
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  13 300     ${OUT_DIR}/prefit-Buenrostro2018-k=13
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  14 300     ${OUT_DIR}/prefit-Buenrostro2018-k=14
sbatch --mem=20G ${SCRIPT_PREFIT} ${DAT_DIR}/Buenrostro_2018_counts.RData  15 300     ${OUT_DIR}/prefit-Buenrostro2018-k=15


# Fit factorizations to the Buenrostro_2018 counts data
#                              data                                    prefitfile                            k  method numiter ex  outfile
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=2  2  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=2
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=3  3  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=3
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=4  4  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=4
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=5  5  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=5
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=6  6  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=6
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=7  7  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=7
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=8  8  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=8
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=9  9  scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=9
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=10 10 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=10
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=11 11 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=11
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=12 12 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=12
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=13 13 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=13
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=14 14 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=14
sbatch --mem=20G ${SCRIPT_FIT} ${DAT_DIR}/Buenrostro_2018_counts.RData ${OUT_DIR}/prefit-Buenrostro2018-k=15 15 scd    600     yes ${OUT_DIR}/fit-Buenrostro2018-scd-ex-k=15

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/unbinarized/
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Buenrostro2018
