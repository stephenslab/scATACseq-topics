#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for Lareau 2019 single-cell ATAC-seq data sets, and
# for different choices of the model parameters and optimization
# settings.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch

# Lareau_2019_bonemarrow data
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/bone_marrow/processed_data/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/bone_marrow/output/
cd /project2/mstephens/kevinluo/scATACseq-topics/data/log
mkdir -p ${OUT_DIR}
NUMITER=10

# "Pre-fit" factorizations to the Lareau_2019_bonemarrow data.
#                       data                                     k numiter    outfile
sbatch ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData  2 ${NUMITER} ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2

# Fit factorizations to Lareau_2019_bonemarrow data, with and without extrapolation.
#                    data                                    prefitfile                                  k method numiter     ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 em     ${NUMITER}  no ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 ccd    ${NUMITER}  no ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 scd    ${NUMITER}  no ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 em     ${NUMITER} yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 ccd    ${NUMITER} yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_bonemarrow.RData ${OUT_DIR}/prefit-Lareau2019_bonemarrow-k=2 2 scd    ${NUMITER} yes ${OUT_DIR}/fit-Lareau2019_bonemarrow-scd-ex-k=2


# Lareau_2019_mousebrain data
DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/mouse_brain/processed_data/
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Lareau_2019/mouse_brain/output/
cd /project2/mstephens/kevinluo/scATACseq-topics/data/log
mkdir -p ${OUT_DIR}
NUMITER=10

# "Pre-fit" factorizations to the Lareau_2019_mousebrain data.
#                       data                                     k numiter    outfile
sbatch ${SCRIPT_PREFIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData  2 ${NUMITER} ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2

# Fit factorizations to Lareau_2019_mousebrain data, with and without extrapolation.
#                    data                                    prefitfile                                  k method numiter     ex outfile
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 em     ${NUMITER}  no fit-Lareau2019_mousebrain-em-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 ccd    ${NUMITER}  no fit-Lareau2019_mousebrain-ccd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 scd    ${NUMITER}  no fit-Lareau2019_mousebrain-scd-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 em     ${NUMITER} yes fit-Lareau2019_mousebrain-em-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 ccd    ${NUMITER} yes fit-Lareau2019_mousebrain-ccd-ex-k=2
sbatch ${SCRIPT_FIT} ${DAT_DIR}/Lareau_2019_mousebrain.RData ${OUT_DIR}/prefit-Lareau2019_mousebrain-k=2 2 scd    ${NUMITER} yes fit-Lareau2019_mousebrain-scd-ex-k=2


