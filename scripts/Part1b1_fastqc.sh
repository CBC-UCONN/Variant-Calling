#!/bin/bash 
#SBATCH --job-name=fastqc
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
module load fastqc

# create output directory
OUTDIR=../results/fastqc
mkdir -p $OUTDIR

# run fastqc. "*fq" tells it to run on the illumina fastq files in directory "../rawdata/"
fastqc -t 6 -o $OUTDIR ../rawdata/*1.fq.gz
fastqc -t 6 -o $OUTDIR ../rawdata/*2.fq.gz