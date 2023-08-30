#!/bin/bash 
#SBATCH --job-name=freebayes
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

# load required software

module load bedtools/2.29.0
module load bamtools/2.5.1
module load htslib/1.12
module load freebayes/1.3.4


# directories/files

INDIR=../../results/03_AlignmentAndCoverage/bwa_align/

OUTDIR=../../results/04_VariantCalling_shortread/variants_freebayes
mkdir -p $OUTDIR 

# make a list of bam files
ls $INDIR/*.bam >$INDIR/bam_list

# set a variable for the reference genome location
GEN=../../genome/GRCh38_latest_genomic.fna

# optionally define a targets file for variant calling in specific regions. 
TARGETS=../../results/03_AlignmentAndCoverage/short_read_coverage/targets.bed.gz

# call freebayes
	# coverage limits could also be set in freebayes using the following flags, e.g.:
	# --min-coverage 110 \
	# --skip-coverage 330 \

freebayes \
-f $GEN \
--bam-list $INDIR/bam_list \
-m 30 \
-q 20 \
-t <(zcat $TARGETS) | \
bgzip -c >$OUTDIR/ashtrio_fb.vcf.gz


tabix -p vcf $OUTDIR/ashtrio_fb.vcf.gz

date

