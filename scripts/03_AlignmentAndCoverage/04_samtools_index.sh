#!/bin/bash
#SBATCH --job-name=samtools_index
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date


module load samtools/1.16.1

INDEXDIR=../../results/03_AlignmentAndCoverage/samtools_index
mkdir -p $INDEXDIR


GENOME=../../genome/GRCh38_latest_genomic.fna

samtools faidx $GENOME -o $INDEXDIR/GRCh38_latest_genomic.fna.fai
