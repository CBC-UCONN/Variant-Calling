#!/bin/bash 
#SBATCH --job-name=gatk_HC
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
module load htslib
module load bedtools

INDIR=../align_pipe

OUTDIR=../variants_gatk
mkdir -p $OUTDIR

# make a list of bam files
ls $INDIR/*bam >$OUTDIR/bam.list

TARGETS=../coverage_stats/targets.bed.gz

# set a variable for the reference genome location
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

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

