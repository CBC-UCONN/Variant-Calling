#!/bin/bash 
#SBATCH --job-name=bcf_tabix
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

module load htslib/1.16

# set input directory
INDIR=../../results/04_VariantCalling_shortread/variants_bcftools

tabix -p vcf $INDIR/ashtrio.vcf.gz

date
