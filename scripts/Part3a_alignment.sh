#!/bin/bash 
#SBATCH --job-name=align_pipe
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
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
module load samtools
module load samblaster
module load bwa/0.7.17

# raw data directory
INDIR=../rawdata

# specify and create output directory
OUTDIR=../align_pipe
mkdir -p $OUTDIR

# set a variable 'GEN' that gives the location and base name of the reference genome:
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38

# execute the pipe for the son:
bwa mem -t 7 -R '@RG\tID:son\tSM:son' $GEN $INDIR/son.1.fq $INDIR/son.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$USER - >$OUTDIR/son.bam
date

# execute the pipe for the mom:
bwa mem -t 7 -R '@RG\tID:mom\tSM:mom' $GEN $INDIR/mom.1.fq $INDIR/mom.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$USER - >$OUTDIR/mom.bam
date

# execute the pipe for the dad:
bwa mem -t 7 -R '@RG\tID:dad\tSM:dad' $GEN $INDIR/dad.1.fq $INDIR/dad.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$USER - >$OUTDIR/dad.bam
date


