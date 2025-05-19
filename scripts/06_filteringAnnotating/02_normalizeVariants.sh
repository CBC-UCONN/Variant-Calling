#!/bin/bash 
#SBATCH --job-name=normalizeVariants
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
module load vcflib/1.0.0-rc1

# vcf files
FREEBAYES=../../results/05_variantCalling/freebayes/freebayes_filtered.vcf.gz
FREEBAYESOUT=../../results/05_variantCalling/freebayes/freebayes_norm.vcf.gz
FREEBAYESAP=../../results/05_variantCalling/freebayes/freebayes_normAP.vcf.gz

GATK=../../results/05_variantCalling/gatk/gatk_filtered.vcf.gz
GATKOUT=../../results/05_variantCalling/gatk/gatk_norm.vcf.gz
GATKOUTAP=../../results/05_variantCalling/gatk/gatk_normAP.vcf.gz

BCFTOOLS=../../results/05_variantCalling/bcftools/bcftools_filtered.vcf.gz
BCFTOOLSOUT=../../results/05_variantCalling/bcftools/bcftools_norm.vcf.gz
BCFTOOLSOUTAP=../../results/05_variantCalling/bcftools/bcftools_normAP.vcf.gz

# genome
GENOME=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

# add vcfstreamsort because bcftools has at least one out of order record after vcfallelicprimitives

# freebayes
bcftools view ${FREEBAYES} | 
  bcftools norm -f ${GENOME} -O z -o ${FREEBAYESOUT}

bcftools view ${FREEBAYES} | 
  bcftools norm -f ${GENOME} |
  vcfallelicprimitives | vcfstreamsort |
  bgzip >${FREEBAYESAP}
tabix -p vcf ${FREEBAYESAP}

# gatk
bcftools view ${GATK} | 
  bcftools norm -f ${GENOME} -O z -o ${GATKOUT}

bcftools view ${GATK} | 
  bcftools norm -f ${GENOME} |
  vcfallelicprimitives | vcfstreamsort |
  bgzip >${GATKOUTAP}
tabix -p vcf ${GATKOUTAP}

# bcftools
bcftools view ${BCFTOOLS} | 
  bcftools norm -f ${GENOME} -O z -o ${BCFTOOLSOUT}

bcftools view ${BCFTOOLS} | 
  bcftools norm -f ${GENOME} |
  vcfallelicprimitives | vcfstreamsort |
  bgzip >${BCFTOOLSOUTAP}
tabix -p vcf ${BCFTOOLSOUTAP}

