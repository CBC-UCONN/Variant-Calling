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

OUTDIR=../results/variants_gatk

# set a variable for the reference genome location
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk GenotypeGVCFs \
    -R $GEN \
    -V gendb://../results/variants_genomicsdb \
    -O $OUTDIR/ashtrio.vcf 

bgzip $OUTDIR/ashtrio.vcf 
tabix -p vcf $OUTDIR/ashtrio.vcf.gz

date