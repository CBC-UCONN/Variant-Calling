#!/bin/bash 
#SBATCH --job-name=bwa_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=20G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load bwa/0.7.17

# set directories

INDEXDIR=../../results/03_AlignmentAndCoverage/bwa_index
mkdir -p $INDEXDIR

GENOME=../../genome/GRCh38_latest_genomic.fna


# Create index for genome
bwa index \
   -p $INDEXDIR/GRCh38 \
   $GENOME


