#!/bin/bash 
#SBATCH --job-name=index_bams
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

module load samtools/1.12

# "*mkdup.bam" will refer to each of the 
for file in ../results/align_stepwise/*mkdup.bam
	do samtools index $file
done

date