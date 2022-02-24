#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for each tissue in the Cusanovich 2018 single-cell ATAC-seq data sets.
SCRIPT_PREFIT=~/projects/scATACseq-topics/scripts/prefit_poisson_nmf.sbatch
SCRIPT_FIT=~/projects/scATACseq-topics/scripts/fit_poisson_nmf.sbatch

# ============ heart, kidney and lung data ============

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_HeartKidneyLung.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_HeartKidneyLung
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the Cusanovich_2018 data.

#                                                    data        k  numiter outfile
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                           k method numiter ex  outfile
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018

# ============ bone marrow data ============

mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_BoneMarrow.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_BoneMarrow
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the bone marrow data.
#                                 data        k  numiter outfile
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to bone marrow data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                           k method numiter ex  outfile
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018


# ============ whole brain ============
mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_WholeBrain.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_WholeBrain
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the whole brain data.
#                                                    data        k  numiter outfile
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to whole brain data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                           k method numiter ex  outfile
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018

# ============ PreFrontalCortex ============
mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_PreFrontalCortex.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_PreFrontalCortex
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the pre-frontal cortex data.
#                                                    data        k  numiter outfile
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to whole brain data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                           k method numiter ex  outfile
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
sbatch --mem=40G --account=pi-xinhe ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018

# ============ Kidney ============
mkdir -p /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_Kidney.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_Kidney
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the whole brain data.
#                                 data        k  numiter outfile
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
sbatch --mem=40G ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanovich2018-k=13

# Fit factorizations to whole brain data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                           k method numiter ex  outfile
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2 2 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3 3 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4 4 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5 5 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6 6 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7 7 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8 8 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9 9 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
sbatch --mem=40G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018

# ============ Endothelial cells ============
cd /project2/mstephens/kevinluo/scATACseq-topics/log/Cusanovich_2018

DAT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data
DATFILE=${DAT_DIR}/Cusanovich_2018_Endothelial.RData
OUT_DIR=/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018_Endothelial
mkdir -p ${OUT_DIR}

# "Pre-fit" factorizations to the endothelial cells data.

#                                 data        k  numiter outfile
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  2  200     ${OUT_DIR}/prefit-Cusanovich2018-k=2
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  3  200     ${OUT_DIR}/prefit-Cusanovich2018-k=3
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  4  200     ${OUT_DIR}/prefit-Cusanovich2018-k=4
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  5  200     ${OUT_DIR}/prefit-Cusanovich2018-k=5
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  6  200     ${OUT_DIR}/prefit-Cusanovich2018-k=6
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  7  200     ${OUT_DIR}/prefit-Cusanovich2018-k=7
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  8  200     ${OUT_DIR}/prefit-Cusanovich2018-k=8
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  9  200     ${OUT_DIR}/prefit-Cusanovich2018-k=9
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  10 200     ${OUT_DIR}/prefit-Cusanovich2018-k=10
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  11 200     ${OUT_DIR}/prefit-Cusanovich2018-k=11
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  12 200     ${OUT_DIR}/prefit-Cusanovich2018-k=12
# sbatch --mem=50G ${SCRIPT_PREFIT} ${DATFILE}  13 200     ${OUT_DIR}/prefit-Cusanotvich2018-k=13

# Fit factorizations to Cusanovich_2018 data, with and without extrapolation.
# Computation took 105612 seconds (29 hrs) for of 250 iterations with k = 7
#                              data       prefitfile                            k  method numiter ex  outfile
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=2  2  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=2
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=3  3  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=3
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=4  4  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=4
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=5  5  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=5
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=6  6  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=6
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=7  7  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=7
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=8  8  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=8
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=9  9  scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=9
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=10 10 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=10
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=11 11 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=11
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=12 12 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=12
# sbatch --mem=50G ${SCRIPT_FIT} ${DATFILE} ${OUT_DIR}/prefit-Cusanovich2018-k=13 13 scd    250     yes ${OUT_DIR}/fit-Cusanovich2018-scd-ex-k=13

# Compile the fitted Poisson non-negative factorizations into a single .RData file.
# Rscript ~/projects/scATACseq-topics/scripts/compile_poisson_nmf_fits.R -o ${OUT_DIR} -d Cusanovich2018
