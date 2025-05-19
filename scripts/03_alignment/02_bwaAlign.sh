#!/bin/bash 
#SBATCH --job-name=align_pipe
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[0-2]

hostname
date

# load required software
module load samtools/1.16.1
module load samblaster/0.1.24
module load bwa-mem2/2.1

#set directories
SAMPDIR=../../results/02_qc/trimmed_fastq

OUTDIR=../../results/03_Alignment/bwa_align
mkdir -p $OUTDIR

INDEX=../../results/03_Alignment/bwa_index/GRCh38

# sample ID list
SAMPLELIST=(son dad mom)

# extract one sample ID
SAMPLE=${SAMPLELIST[$SLURM_ARRAY_TASK_ID]}

# create read group string
RG=$(echo \@RG\\tID:$SAMPLE\\tSM:$SAMPLE)

# execute the alignment pipe:
bwa-mem2 mem -t 7 -R ${RG} ${INDEX} ${SAMPDIR}/${SAMPLE}_trim.1.fq.gz $SAMPDIR/${SAMPLE}_trim.2.fq.gz | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/${SAMPLE}.temp -O BAM >$OUTDIR/${SAMPLE}.bam 

# index alignment file
samtools index ${OUTDIR}/${SAMPLE}.bam

