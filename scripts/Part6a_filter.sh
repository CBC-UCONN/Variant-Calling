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

module load bcftools/1.12
module load htslib/1.12
module load bedtools/2.29.0

OUTDIR=../results/filtered_vcfs
mkdir -p $OUTDIR

TARGETS=../results/coverage_stats/targets.bed

# using a subset (chr20:31400000-34400000) which excludes the most problematic areas

# jointly called variants for illumina data -------------------------

bcftools view -r chr20:31400000-34400000 ../results/variants_freebayes/ashtrio_fb.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/fb_ill_filter.vcf.gz

bcftools view -r chr20:31400000-34400000 ../results/variants_gatk/ashtrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/gatk_ill_filter.vcf.gz
# additionally filter bcftools output to match targets used for gatk and freebayes
bcftools view -r chr20:31400000-34400000 ../results/variants_bcftools/ashtrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bedtools intersect -header -wa -a stdin -b $TARGETS | bgzip -c > $OUTDIR/bcf_ill_filter.vcf.gz

# son only to compare long reads against short reads for a single sample -------------------------

bcftools view -r chr20:31400000-34400000 ../results/variants_clair3_gvcf/son/merge_output.vcf.gz | bedtools intersect -header -wa -a stdin -b $TARGETS | bgzip -c > $OUTDIR/clair3_son_ont_filter.vcf.gz

bcftools view -r chr20:31400000-34400000 ../results/variants_pepper/son/son.vcf.gz | bedtools intersect -header -wa -a stdin -b $TARGETS | bgzip -c > $OUTDIR/pepper_son_ont_filter.vcf.gz

bcftools view -r chr20:31400000-34400000 ../results/variants_gatk/son.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bedtools intersect -header -wa -a stdin -b $TARGETS | bgzip -c > $OUTDIR/gatk_son_ont_filter.vcf.gz

# index vcfs -------------------

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done
