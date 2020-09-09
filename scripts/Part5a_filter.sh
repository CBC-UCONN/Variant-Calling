#!/bin/bash
#SBATCH --job-name=bcf_vcf
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools/1.6
module load htslib

OUTDIR=../filtered_vcfs
mkdir -p $OUTDIR

bcftools filter -s LowQual -e '%QUAL<50' ../variants_freebayes/chinesetrio_fb.vcf.gz | bgzip -c > $OUTDIR/fb_filter.vcf.gz
bcftools filter -s LowQual -e '%QUAL<50' ../variants_gatk/chinesetrio.vcf.gz | bgzip -c > $OUTDIR/gatk_filter.vcf.gz
bcftools filter -s LowQual -e '%QUAL<50' ../variants_bcftools/chinesetrio.vcf.gz | bgzip -c > $OUTDIR/bcf_filter.vcf.gz

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done



