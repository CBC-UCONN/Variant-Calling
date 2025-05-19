#!/bin/bash 
#SBATCH --job-name=genmap_mappable
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --mem=80G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load genmap/1.3.0

# define and/or create input, output directories
INDEXDIR=../../results/04_alignQC/genmapindex

OUTDIR=../../results/04_alignQC/genmapmappability
mkdir -p ${OUTDIR}

PREFIX=mappable

# run mappability algorithm
genmap map -K 30 -E 2 -I ${INDEXDIR} -O ${OUTDIR}/${PREFIX} -t -w -bg
