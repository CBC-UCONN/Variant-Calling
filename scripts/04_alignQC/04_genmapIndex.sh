#!/bin/bash 
#SBATCH --job-name=genmap_index
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
GENOME=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

OUTDIR=../../results/04_alignQC/
mkdir -p ${OUTDIR}

INDEXDIR=genmapindex

genmap index -F ${GENOME} -I ${OUTDIR}/${INDEXDIR}
