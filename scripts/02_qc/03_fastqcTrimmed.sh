#!/bin/bash 
#SBATCH --job-name=fastqc_trimmed
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

# input/output directory
INDIR=../../results/02_qc/trimmed_fastq
OUTDIR=../../results/02_qc/fastqc_trimmed
mkdir -p $OUTDIR

# run fastqc. "*fq" tells it to run on the illumina fastq files in directory "data/"
fastqc -t 6 -o $OUTDIR ${INDIR}/*trim.1.fq.gz
fastqc -t 6 -o $OUTDIR ${INDIR}/*trim.2.fq.gz

module load MultiQC/1.9

# run multiqc on fastqc output
multiqc -f -o $OUTDIR/multiqc $OUTDIR