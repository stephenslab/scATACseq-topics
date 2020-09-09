#!/bin/bash

#SBATCH --job-name=fitNMF
#SBATCH --output=fitNMF_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=36:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs fit_poisson_nmf.R. See fit_all_models.sh for examples
# illustrating how to use this script.

# Get the command-line arguments.
DATFILE=${1}
PREFITFILE=${2}
K=${3}
METHOD=${4}
NUMITER=${5}
if [[ ${6} == "yes" ]]
then
  EXTRAPOLATE="--extrapolate"
else
  EXTRAPOLATE=""
fi
OUTFILE=${7}
echo datfile=${DATFILE}
echo prefitfile=${PREFITFILE}
echo k=${K}
echo method=${METHOD}
echo numiter=${NUMITER}
echo extrapolate=${6}
echo outfile=${OUTFILE}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1

# Run the R script.
echo "Running fit_poisson_nmf.R"
export OPENBLAS_NUM_THREADS=1
Rscript ~/projects/scATACseq-topics/scripts/fit_poisson_nmf.R \
  --counts ${DATFILE} --prefit ${PREFITFILE}.rds \
  --method ${METHOD} -n ${NUMITER} ${EXTRAPOLATE} \
  --nc 10 -o ${OUTFILE}.rds > ${OUTFILE}.Rout
