#!/bin/bash 
#SBATCH --job-name=filterVariants
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

module load bcftools/1.16
module load htslib/1.16
module load bedtools/2.29.0

# VCF files
FREEBAYES=../../results/05_variantCalling/freebayes/freebayes.vcf.gz
FREEBAYESOUT=../../results/05_variantCalling/freebayes/freebayes_filtered.vcf.gz

GATK=../../results/05_variantCalling/gatk/gatk.vcf.gz
GATKOUT=../../results/05_variantCalling/gatk/gatk_filtered.vcf.gz

BCFTOOLS=../../results/05_variantCalling/bcftools/bcftools.vcf.gz
BCFTOOLSOUT=../../results/05_variantCalling/bcftools/bcftools_filtered.vcf.gz

# Targets file
COV=../../results/04_alignQC/coverage/coverage_1kb.bed.gz
TARGETS=../../results/04_alignQC/coverage/targets.bed.gz

zcat ${COV} |
awk '$4 > 50 && $4 < 250' |
bedtools merge -i stdin | 
gzip >${TARGETS}

# filtering

# freebayes: use AB/AF/DP/QUAL
bcftools filter \
  -r chr20:29400000-34400000 \
  --exclude 'INFO/DP < 50 | INFO/DP > 250 | AB < .25 | AB > 0.75 | QUAL < 30 | AF < 0.1' \
  ${FREEBAYES} |
bedtools intersect -header -a stdin -b ${TARGETS} |
bgzip >${FREEBAYESOUT}

tabix -p vcf ${FREEBAYESOUT}

# gatk: use DP/QUAL (others not available)
bcftools filter \
  -r chr20:29400000-34400000 \
  --exclude 'INFO/DP < 50 | INFO/DP > 250 | QUAL < 30' \
  ${GATK} |
bedtools intersect -header -a stdin -b ${TARGETS} |
bgzip >${GATKOUT}

tabix -p vcf ${GATKOUT}

# bcftools: use DP/QUAL (others not available)
bcftools filter \
  -r chr20:29400000-34400000 \
  --exclude 'INFO/DP < 50 | INFO/DP > 250 | QUAL < 30' \
  ${BCFTOOLS} |
bedtools intersect -header -a stdin -b ${TARGETS} |
bgzip >${BCFTOOLSOUT}

tabix -p vcf ${BCFTOOLSOUT}