#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################

module load fastp/0.23.2

INDIR=../../data

REPORTDIR=../../results/02_QualityControl/fastp_reports
mkdir -p $REPORTDIR

TRIMDIR=../../results/02_QualityControl/trimmed
mkdir -p $TRIMDIR


for infile in $INDIR/*1.fq.gz


do


base=$(basename ${infile} .1.fq.gz)


fastp \
    --in1 ${infile} \
    --in2 $INDIR/${base}.2.fq.gz \
    --out1 $TRIMDIR/${base}.trim.1.fq \
    --out2 $TRIMDIR/${base}.trim.2.fq \
    --json $REPORTDIR/${base}_fastp.json \
    --html $REPORTDIR/${base}_fastp.html

done

