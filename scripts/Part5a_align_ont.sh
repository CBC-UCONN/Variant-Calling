#!/bin/bash
#SBATCH --job-name=align_ont
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


# load software

module load minimap2/2.18
module load samtools/1.12

# input/output files, directories
INDIR=../rawdata

OUTDIR=../results/align_ont/
mkdir -p $OUTDIR

NPROC=$(nproc)
GENOME=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta


# run minimap

# son 
SAM=son
minimap2 -c --MD -ax map-ont -t 10 $GENOME $INDIR/$SAM.ont.fq.gz | \
samtools sort -@ 5 -T $OUTDIR/$SAM.temp -O BAM \
>$OUTDIR/$SAM.bam

samtools index $OUTDIR/$SAM.bam

# mother
SAM=mom
minimap2 -c --MD -ax map-ont -t 10 $GENOME $INDIR/$SAM.ont.fq.gz | \
samtools sort -@ 5 -T $OUTDIR/$SAM.temp -O BAM \
>$OUTDIR/$SAM.bam

samtools index $OUTDIR/$SAM.bam

# father
SAM=dad
minimap2 -c --MD -ax map-ont -t 10 $GENOME $INDIR/$SAM.ont.fq.gz | \
samtools sort -@ 5 -T $OUTDIR/$SAM.temp -O BAM \
>$OUTDIR/$SAM.bam

samtools index $OUTDIR/$SAM.bam
