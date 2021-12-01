#!/bin/bash 
#SBATCH --job-name=compress
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
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
module load samtools/1.12

# set directory for input and output files
DIR=../results/align_stepwise

# son
samtools view -bhS $DIR/son.sam >$DIR/son.bam
# mom
samtools view -bhS $DIR/mom.sam >$DIR/mom.bam
# dad
samtools view -bhS $DIR/dad.sam >$DIR/dad.bam

date