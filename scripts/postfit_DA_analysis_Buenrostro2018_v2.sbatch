#!/bin/bash

#SBATCH --job-name=DA_analysis
#SBATCH --output=DA_analysis_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G

# For reproducibility, I add version numbers to the "module load" calls.
module load R/4.0.4

# Run the R script.
echo "Running DA_analysis_Buenrostro2018.R"
Rscript ~/projects/scATACseq-topics/scripts/DA_analysis_Buenrostro2018.R
