#!/bin/bash 
#SBATCH --job-name=bcf_mpileup
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools

INDIR=../align_stepwise
# make output directory if it doesn't exist. 
OUTDIR=../variants_bcftools
mkdir -p $OUTDIR


# make a list of bam files
ls $INDIR/*mkdup.bam >$INDIR/list.bam

# set reference genome location
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

bcftools mpileup \
	-f $GEN \
	-b $INDIR/list.bam \
	-q 20 \
	-Q 30 \
	-r chr20:29400000-34400000 >$OUTDIR/chinesetrio.pileup
date

