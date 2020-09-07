#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load bwa/0.7.17
module load samtools

# set input and output directories
INDIR=../rawdata
OUTDIR=../align_stepwise

mkdir -p $OUTDIR


# current location of indexed HG38
# may need to be changed. 
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38

# note that read group info is added during alignment. 

# each line aligns one family member's sequences
# son
bwa mem -t 12 -R '@RG\tID:son\tSM:son' $GEN $INDIR/son.1.fq $INDIR/son.2.fq -o $OUTDIR/son.sam 
date
# mom
bwa mem -t 12 -R '@RG\tID:mom\tSM:mom' $GEN $INDIR/mom.1.fq $INDIR/mom.2.fq -o $OUTDIR/mom.sam 
date
# dad
bwa mem -t 12 -R '@RG\tID:dad\tSM:dad' $GEN $INDIR/dad.1.fq $INDIR/dad.2.fq -o $OUTDIR/dad.sam 
date

