#!/bin/bash
#SBATCH --job-name=filter_vcfs
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
module load bedtools/2.29.0

OUTDIR=../filtered_vcfs
mkdir -p $OUTDIR

TARGETS=../coverage_stats/targets.bed

# using a subset (chr20:31500000-34400000) for further analysis b/c mom/dad and son bam files on the giab ftp site are on different reference genomes, so there is a bunch of missing data. 

bcftools view -r chr20:31500000-34400000 ../variants_freebayes/chinesetrio_fb.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/fb_filter.vcf.gz
bcftools view -r chr20:31500000-34400000 ../variants_gatk/chinesetrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/gatk_filter.vcf.gz
# additionally filter bcftools output to match targets used for gatk and freebayes
bcftools view -r chr20:31500000-34400000 ../variants_bcftools/chinesetrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bedtools intersect -header -wa -a stdin -b $TARGETS | bgzip -c > $OUTDIR/bcf_filter.vcf.gz

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done



