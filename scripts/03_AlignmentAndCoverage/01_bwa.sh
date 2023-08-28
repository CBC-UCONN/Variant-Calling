#!/bin/bash 
#SBATCH --job-name=align_pipe
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=20G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load samtools/1.16.1
module load samblaster/0.1.24
module load bwa/0.7.17

#set directories
SAMPDIR=../../results/02_QualityControl/trimmed


OUTDIR=../../results/03_AlignmentAndCoverage/bwa_align
mkdir -p $OUTDIR

INDEXDIR=../../results/03_AlignmentAndCoverage/bwa_index
mkdir -p $INDEXDIR

GENOME=../../genome/GRCh38_latest_genomic.fna


###First we need to create index for genome
#bwa index \
 #   -p $INDEXDIR/GRCh38 \
  #  $GENOME

INDEX=../../results/03_AlignmentAndCoverage/bwa_index/GRCh38

# execute the pipe for the son:
bwa mem -t 7 -R '@RG\tID:son\tSM:son' $INDEX $SAMPDIR/son.trim.1.fq $SAMPDIR/son.trim.2.fq | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/son.temp -O BAM >$OUTDIR/son.bam 

samtools stats ${OUTDIR}/son.bam >${OUTDIR}/son.stats


# execute the pipe for the mom:
bwa mem -t 7 -R '@RG\tID:mom\tSM:mom' $INDEX $SAMPDIR/mom.trim.1.fq $SAMPDIR/mom.trim.2.fq | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/mom.temp -O BAM >$OUTDIR/mom.bam 

samtools stats ${OUTDIR}/mom.bam >${OUTDIR}/mom.stats


# execute the pipe for the dad:
bwa mem -t 7 -R '@RG\tID:dad\tSM:dad' $INDEX $SAMPDIR/dad.trim.1.fq $SAMPDIR/dad.trim.2.fq | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/dad.temp -O BAM >$OUTDIR/dad.bam 

samtools stats ${OUTDIR}/dad.bam >${OUTDIR}/dad.stats
