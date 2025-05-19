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
OUTDIR=../../data/genomeWide
mkdir -p $OUTDIR

# son -----------------------------
SONR1='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_L002_R1_003.fastq.gz'
SONR2='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_L002_R2_003.fastq.gz'

wget -P ${OUTDIR} ${SONR1}
wget -P ${OUTDIR} ${SONR2}

# mom -----------------------------
MOMR1='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/reads/D3_S3_L001_R1_004.fastq.gz'
MOMR2='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/reads/D3_S3_L001_R2_004.fastq.gz'

wget -P ${OUTDIR} ${MOMR1}
wget -P ${OUTDIR} ${MOMR2}

# dad -----------------------------
DADR1='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/reads/D2_S2_L002_R1_001.fastq.gz'
DADR2='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/reads/D2_S2_L002_R2_001.fastq.gz'

wget -P ${OUTDIR} ${DADR1}
wget -P ${OUTDIR} ${DADR2}