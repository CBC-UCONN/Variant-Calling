#!/bin/bash 
#SBATCH --job-name=alignPipe
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
DATA=../../data/genomeWide/

OUTDIR=../../results/07_parallelizing/bwaAlign
mkdir -p ${OUTDIR}

INDEX=../../results/03_Alignment/bwa_index/GRCh38

# get sample ID, R1, R2 from samples.csv
SN=$( expr ${SLURM_ARRAY_TASK_ID} + 1 )
SAMPLE=$(cut -d "," -f 1 samples.csv | sed -n ${SN}p)
R1=$(cut -d "," -f 2 samples.csv | sed -n ${SN}p)
R2=$(cut -d "," -f 3 samples.csv | sed -n ${SN}p)

# create read group string
RG=$(echo \@RG\\tID:$SAMPLE\\tSM:$SAMPLE)

# execute the alignment pipe:
bwa-mem2 mem -t 7 -R ${RG} ${INDEX} ${DATA}/${R1} $DATA/${R2} | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/${SAMPLE}.temp -O BAM >$OUTDIR/${SAMPLE}.bam 

# index alignment file
samtools index ${OUTDIR}/${SAMPLE}.bam

