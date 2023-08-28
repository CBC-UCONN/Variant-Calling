#!/bin/bash 
#SBATCH --job-name=index_bwa
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load samtools/1.16.1

for file in ../../results/03_AlignmentAndCoverage/bwa_align/*bam
	do samtools index $file
done
