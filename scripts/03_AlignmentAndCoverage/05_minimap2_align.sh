#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software

module load minimap2/2.24
module load samtools/1.16.1

# input/output files, directories
ONT=../../data

OUTDIR=../../results/03_AlignmentAndCoverage/minimap2
    mkdir -p ${OUTDIR}

GENOME=../../genome/GRCh38_latest_genomic.fna


#run minimap

for infile in $ONT/*ont.fq.gz


do

OUTROOT=$(basename ${infile} .ont.fq.gz)


# run minimap
minimap2 -c --MD -ax map-ont -t 15 ${GENOME} ${infile} | \
samtools sort -@ 5 -T ${OUTDIR}/${OUTROOT}.temp -O BAM >${OUTDIR}/${OUTROOT}.bam

samtools index ${OUTDIR}/${OUTROOT}.bam

samtools stats ${OUTDIR}/${OUTROOT}.bam > ${OUTDIR}/${OUTROOT}.stats


done
