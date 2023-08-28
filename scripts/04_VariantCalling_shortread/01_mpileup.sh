#!/bin/bash 
#SBATCH --job-name=bcf_mpileup
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

module load bcftools/1.16

INDIR=../../results/03_AlignmentAndCoverage/bwa_align
# make output directory if it doesn't exist. 
OUTDIR=../../results/04_VariantCalling_shortread/variants_bcftools
mkdir -p $OUTDIR


# make a list of bam files
ls $INDIR/*.bam >$INDIR/bam_list

# set reference genome location
GEN=../../genome/GRCh38_latest_genomic.fna


bcftools mpileup -f $GEN -b $INDIR/bam_list -q 20 -Q 30 | bcftools call -m -v -Oz -o $OUTDIR/ashtrio.vcf.gz 

date
