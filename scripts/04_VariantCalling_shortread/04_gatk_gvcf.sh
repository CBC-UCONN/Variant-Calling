#!/bin/bash 
#SBATCH --job-name=gatk_gvcf
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
module load picard/2.23.9

INDIR=../../results/03_AlignmentAndCoverage/bwa_align/

OUTDIR=../../results/04_VariantCalling_shortread/variants_gatk
mkdir -p $OUTDIR

# make a list of bam files
ls $INDIR/*.bam >$INDIR/bam_list

# optionally define a targets file for variant calling in specific regions. 
TARGETS=../../results/03_AlignmentAndCoverage/short_read_coverage/targets.bed.gz

# set a variable for the reference genome location
GEN=../../genome/GRCh38_latest_genomic.fna

# create required .dict file
java -jar $PICARD CreateSequenceDictionary R=$GEN

# run haplotype caller on each sample
gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/son.bam \
     -ERC GVCF \
     -L $TARGETS \
     --output $OUTDIR/son.g.vcf

date 

gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/mom.bam \
     -ERC GVCF \
     -L $TARGETS \
     --output $OUTDIR/mom.g.vcf

date 

gatk HaplotypeCaller \
     -R $GEN \
     -I $INDIR/dad.bam \
     -ERC GVCF \
     -L $TARGETS \
     --output $OUTDIR/dad.g.vcf

date

