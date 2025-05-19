#!/bin/bash 
#SBATCH --job-name=genotypeGVCFs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load GATK/4.0
module load htslib/1.16

OUTDIR=../../results/05_variantCalling/gatk


# set a variable for the reference genome location
GEN=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

gatk GenotypeGVCFs \
    -R ${GEN} \
    -V gendb://../../results/05_variantCalling/gatk/db \
    -O ${OUTDIR}/gatk.vcf 

bgzip ${OUTDIR}/gatk.vcf 
tabix -p vcf ${OUTDIR}/gatk.vcf.gz

