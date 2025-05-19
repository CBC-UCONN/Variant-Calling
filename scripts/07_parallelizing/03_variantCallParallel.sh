#!/bin/bash
#SBATCH --job-name=bcftoolsParallelized
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=100G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#Load software
module load bcftools/1.20
module load htslib/1.20
module load samtools/1.20
module load parallel/20180122
module load vcflib/1.0.0-rc1
module load bedtools/2.29.0

#Specify input/output directories
INDIR=../../results/07_parallelizing/bwaAlign
OUTDIR=../../results/07_parallelizing/bcftoolsVariants
	mkdir -p ${OUTDIR}

#Make a list of bam files
ls ${INDIR}/*.bam >${INDIR}/bam_list.txt

#Set reference genome location
GEN=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

# fai index genome
samtools faidx ${GEN}

#Make regions file
regionsfile=../../genome/100kbWin.bed
bedtools makewindows -g ${GEN}.fai -w 100000 | awk '{start=$2+1}{print $1":"start"-"$3}' >${regionsfile}

#Call variants in parallel
(cat "$regionsfile" | parallel -k -j 20 \
 "bcftools mpileup -f ${GEN} -b ${INDIR}/bam_list.txt -q 25 -Q 25 -r {} | bcftools call -m -v"
) | 
vcffirstheader |
vcfstreamsort -w 1000 |
vcfuniq |
bgzip >${OUTDIR}/bcftools.vcf.gz

#Index vcf
tabix -p vcf ${OUTDIR}/bcftools.vcf.gz