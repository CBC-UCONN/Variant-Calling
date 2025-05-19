#!/bin/bash 
#SBATCH --job-name=freebayes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



hostname
date

# load required software
module load freebayes/1.3.4
module load htslib/1.16

# directories/files

INDIR=../../results/03_Alignment/bwa_align/

OUTDIR=../../results/05_variantCalling/freebayes
mkdir -p ${OUTDIR} 

# make a list of bam files
ls ${INDIR}/*.bam >${INDIR}/bam_list.txt

# set a variable for the reference genome location
GEN=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

# run freebayes
freebayes \
-f $GEN \
--bam-list ${INDIR}/bam_list.txt |
bgzip -c >${OUTDIR}/freebayes.vcf.gz

tabix -p vcf ${OUTDIR}/freebayes.vcf.gz

date