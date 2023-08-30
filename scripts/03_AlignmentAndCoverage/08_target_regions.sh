#!/bin/bash 
#SBATCH --job-name=targets
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



hostname
date

# load required software

module load bedtools/2.29.0

# define and/or create input, output directories

OUTDIR=../../results/03_AlignmentAndCoverage/short_read_coverage/
mkdir -p $OUTDIR


# select and merge outlier windows (after deciding what is an outlier by looking at the distribution in R)
zcat $OUTDIR/coverage_1kb.bed.gz | awk '$4 < 120 || $4 > 250' | bedtools merge | bgzip >$OUTDIR/coverage_outliers.bed.gz 
tabix -p bed $OUTDIR/coverage_outliers.bed.gz

# select and merge target windows (inverse of "outlier windows" above)
zcat $OUTDIR/coverage_1kb.bed.gz | awk '$4 >= 120 && $4 <= 250' | bedtools merge | bgzip >$OUTDIR/targets.bed.gz
tabix -p bed $OUTDIR/targets.bed.gz
