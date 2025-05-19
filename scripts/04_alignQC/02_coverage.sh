#!/bin/bash 
#SBATCH --job-name=coverage
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software

module load bedtools/2.29.0
module load bamtools/2.5.1
module load samtools/1.16.1

# define and/or create input, output directories

INDIR=../../results/03_Alignment/bwa_align

OUTDIR=../../results/04_alignQC/coverage
mkdir -p $OUTDIR

# create faidx genome index file
GENOME=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
FAI=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.fai
samtools faidx ${GENOME}

# make a "genome" file, required by bedtools makewindows command, set variable for location
GFILE=${OUTDIR}/GRCh38.genome
cut -f 1-2 $FAI > $GFILE

# make 1kb window bed file, set variable for location
WIN1KB=${OUTDIR}/GRCh38_1kb.bed
bedtools makewindows -g ${GFILE} -w 1000 >${WIN1KB}

# make a list of bam files
find ${INDIR} -name "*bam" >${OUTDIR}/bam.list

# calculate per-base coverage as well	
bamtools merge -list ${OUTDIR}/bam.list | \
bamtools filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
samtools depth -a /dev/stdin | \
awk '{OFS="\t"}{print $1,$2-1,$2,$3}' | \
bedtools map \
-a ${WIN1KB} \
-b stdin \
-c 4 -o mean,median,count \
-g ${GFILE} | \
bgzip >${OUTDIR}/coverage_1kb.bed.gz