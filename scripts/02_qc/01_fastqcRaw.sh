#!/bin/bash 
#SBATCH --job-name=fastqc_raw
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load fastqc/0.11.7

# create output directory
OUTDIR=../../results/02_qc/fastqc_raw
mkdir -p $OUTDIR

# run fastqc. "*fq" tells it to run on the illumina fastq files in directory "data/"
fastqc -t 6 -o $OUTDIR ../../data/*1.fq.gz
fastqc -t 6 -o $OUTDIR ../../data/*2.fq.gz

module load MultiQC/1.9

# run on fastqc output
multiqc -f -o $OUTDIR/multiqc $OUTDIR