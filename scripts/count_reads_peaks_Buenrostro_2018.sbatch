#!/bin/bash

#SBATCH --job-name=count_reads_peaks
#SBATCH --output=count_reads_peaks_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --mem=20G

# This script counts scATAC-seq reads for ATAC-seq peaks using peak and bam files using `bedtools coverage`
# adapted from https://github.com/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Control/count_reads_peaks.sh

dirinput=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/data/input_Chen_2019
peakfile=$dirinput/combined.sorted.merged.bed
dirout=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/data/count_reads_peaks_output
mkdir -p $dirout

dirBAM=/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/data/input_Chen_2019/sc-bams_nodup/
filecount=0
bamfiles=(`ls $dirBAM/*.bam`)

for bamfile in "${bamfiles[@]}"; do
  echo "input bamfile: $bamfile"
  echo "output: $dirout/$(basename $bamfile).peaks.txt"
  bedtools coverage -a $peakfile -b $bamfile > $dirout/$(basename $bamfile).peaks.txt
  filecount=$((filecount+1))
done

echo "$filecount files counted."
