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


# make a directory for results if it doesn't exist
OUTDIR=../results/variants_freebayes
mkdir -p $OUTDIR 

# make a list of bam files
find ../results/align_pipe/ -name "*bam" >$OUTDIR/bam.list

# set a variable for the reference genome location
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

OUTLIERWINDOWS=../results/coverage_stats/coverage_outliers.bed.gz
TARGETS=../results/coverage_stats/targets.bed

# call freebayes
	# coverage limits could also be defined by looking at the distribution of per base coverage and setting thresholds, e.g.:
	# --min-coverage 110 \
	# --skip-coverage 330 \

freebayes \
-f $GEN \
--bam-list $OUTDIR/bam.list \
-m 30 \
-q 20 \
-t $TARGETS | \
bgzip -c >$OUTDIR/ashtrio_fb.vcf.gz


tabix -p vcf $OUTDIR/ashtrio_fb.vcf.gz

date
