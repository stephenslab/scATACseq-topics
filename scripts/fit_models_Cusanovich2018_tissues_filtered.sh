#!/bin/bash

# The shell commands below will submit slurm jobs to perform the
# Poisson NMF model fitting for each tissue in the Cusanovich 2018 single-cell ATAC-seq data sets.
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf_Cusanovich2018_tissue.sbatch

# ============ PreFrontalCortex ============
mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/tissues/PreFrontalCortex
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/tissues/PreFrontalCortex

#                              tissue            k
sbatch --partition=mstephens --mem=50G ${SCRIPT_FIT} PreFrontalCortex  7
sbatch --mem=50G ${SCRIPT_FIT} PreFrontalCortex  6
sbatch --mem=50G ${SCRIPT_FIT} PreFrontalCortex  8
sbatch --mem=50G ${SCRIPT_FIT} PreFrontalCortex  9
sbatch --mem=50G ${SCRIPT_FIT} PreFrontalCortex  10

# ============ Kidney ============
mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/tissues/Kidney
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018/tissues/Kidney

#                                tissue  k
sbatch --partition=mstephens --mem=40G ${SCRIPT_FIT} Kidney  9
sbatch --mem=40G ${SCRIPT_FIT} Kidney  6
sbatch --mem=40G ${SCRIPT_FIT} Kidney  7
sbatch --mem=40G ${SCRIPT_FIT} Kidney  8
sbatch --mem=40G ${SCRIPT_FIT} Kidney  10
