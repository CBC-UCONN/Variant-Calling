#!/bin/bash 
#SBATCH --job-name=gatk_single_samples
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# make sure partition is specified as `xeon` to prevent slowdowns on amd processors. 

# load required software

module load GATK/4.0
module load htslib/1.12
module load bedtools/2.29.0

INDIR=../../results/03_AlignmentAndCoverage/bwa_align/

OUTDIR=../../results/04_VariantCalling_shortread/variants_gatk_singlesample
mkdir -p $OUTDIR

# make a list of bam files
ls $INDIR/*bam >$OUTDIR/bam.list

TARGETS=../../results/03_AlignmentAndCoverage/short_read_coverage/targets.bed.gz

# set a variable for the reference genome location
GEN=../../genome/GRCh38_latest_genomic.fna

gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/son.bam \
     -L $TARGETS \
     --output $OUTDIR/son.vcf

bgzip $OUTDIR/son.vcf
tabix -p vcf $OUTDIR/son.vcf.gz

gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/mom.bam \
     -L $TARGETS \
     --output $OUTDIR/mom.vcf

bgzip $OUTDIR/mom.vcf 
tabix -p vcf $OUTDIR/mom.vcf.gz 

gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/dad.bam \
     -L $TARGETS \
     --output $OUTDIR/dad.vcf

bgzip $OUTDIR/dad.vcf
tabix -p vcf $OUTDIR/dad.vcf.gz

date

