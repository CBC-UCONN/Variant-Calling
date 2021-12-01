#!/bin/bash
#SBATCH --job-name=trim_sickle
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=10G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname

module load sickle
module load fastqc

# input/output directories
INDIR=../rawdata

OUTDIR=../results/trimmed
mkdir -p $OUTDIR

# trim son
SEQ=son
sickle pe -t sanger \
    -l 100 \
    -f $INDIR/$SEQ.1.fq.gz \
	-r $INDIR/$SEQ.2.fq.gz \
    -o $OUTDIR/$SEQ.trim.1.fq \
    -p $OUTDIR/$SEQ.trim.2.fq \
    -s $OUTDIR/$SEQ.trim.0.fq
# trim mom
SEQ=mom
sickle pe -t sanger \
    -l 100 \
    -f $INDIR/$SEQ.1.fq.gz \
	-r $INDIR/$SEQ.2.fq.gz \
    -o $OUTDIR/$SEQ.trim.1.fq \
    -p $OUTDIR/$SEQ.trim.2.fq \
    -s $OUTDIR/$SEQ.trim.0.fq
# trim dad
SEQ=dad
sickle pe -t sanger \
    -l 100 \
    -f $INDIR/$SEQ.1.fq.gz \
	-r $INDIR/$SEQ.2.fq.gz \
    -o $OUTDIR/$SEQ.trim.1.fq \
    -p $OUTDIR/$SEQ.trim.2.fq \
    -s $OUTDIR/$SEQ.trim.0.fq


# run fastqc on all files matching "../rawdata/*trim*fq"
fastqc -t 6 -o ../results/fastqc/ $OUTDIR/*trim*fq

