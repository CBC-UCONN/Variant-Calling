#!/bin/bash 
#SBATCH --job-name=download_reads
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load samtools/1.12 
module load bedtools/2.29.0

# specify input/output dirs
OUTDIR=../../data/
mkdir -p $OUTDIR

# download a subregion of ashkenazim GIAB trio: chr20:29400000-34400000
# sort reads by name, convert to fastq. 
# data: illumina 2x250bp paired end reads
# there will be lots of errors from bedtools resulting from discordant PE reads. 
	# this is not an issue in this specific tutorial case, but could be worrisome in other contexts. 

# son -----------------------------

# ILLUMINA
SON='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam'
samtools view -uh $SON chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/son.1.fq -fq2 $OUTDIR/son.2.fq

# mom -----------------------------

# ILLUMINA
MOM='ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam'
samtools view -uh $MOM chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/mom.1.fq -fq2 $OUTDIR/mom.2.fq

# dad -----------------------------

# ILLUMINA
DAD='ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam'
samtools view -uh $DAD chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/dad.1.fq -fq2 $OUTDIR/dad.2.fq

# ---------------------------

# get rid of bam indexes that were also downloaded
rm *bam.bai

# gzip fq files
for file in $OUTDIR/*fq
do gzip $file
done