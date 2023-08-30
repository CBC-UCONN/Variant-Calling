#!/bin/bash 
#SBATCH --job-name=gatk_genotype
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
module load htslib/1.12

OUTDIR=../../results/04_VariantCalling_shortread/variants_gatk

# set a variable for the reference genome location
GEN=../../genome/GRCh38_latest_genomic.fna

gatk GenotypeGVCFs \
    -R $GEN \
    -V gendb://../../results/04_VariantCalling_shortread/variants_gatk/db \
    -O $OUTDIR/ashtrio.vcf 

bgzip $OUTDIR/ashtrio.vcf 
tabix -p vcf $OUTDIR/ashtrio.vcf.gz

