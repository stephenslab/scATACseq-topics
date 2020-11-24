#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# cisTopic model fitting for the Buenrostro_2018_binarized_scPeaks data.

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline/chromVAR/

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/cisTopic_Buenrostro_2018_chromVAR_scPeaks
cd /project2/mstephens/kevinluo/scATACseq-topics/log/cisTopic_Buenrostro_2018_chromVAR_scPeaks

# fit cisTopic models with k topics
#------------------------------------------------------------------

OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/cisTopic_Buenrostro_2018_chromVAR_scPeaks/binarized/
mkdir -p ${OUT_DIR}

SCRIPT1=~/projects/scATACseq-topics/scripts/fit_cisTopic.sbatch

#                           data                                                k  n        outfile
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  2  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=2
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  3  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=3
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  4  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=4
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  5  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=5
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  6  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=6
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  7  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=7
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  8  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=8
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  9  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=9
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  10 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=10
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  11 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=11
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  12 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=12
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  13 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=13
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  14 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=14
sbatch --mem=20G ${SCRIPT1} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  15 1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized-k=15


# fit cisTopic models with different numbers of topics
#------------------------------------------------------------------

OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/cisTopic_Buenrostro_2018_chromVAR_scPeaks/binarized/
mkdir -p ${OUT_DIR}

SCRIPT2=~/projects/scATACseq-topics/scripts/fit_cisTopic_topics.sbatch

#                           data                                                n        outfile
sbatch --mem=40G ${SCRIPT2} ${DAT_DIR}/Buenrostro_2018_binarized_scPeaks.RData  1000     ${OUT_DIR}/cisTopic-Buenrostro2018-binarized

