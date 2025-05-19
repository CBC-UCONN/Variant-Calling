#!/bin/bash 
#SBATCH --job-name=bwa_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load required software
module load bwa-mem2/2.1

# set directories
INDEXDIR=../../results/03_Alignment/bwa_index
mkdir -p $INDEXDIR

GENOME=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta


# Create index for genome
    # requires significant memory for indexing
bwa-mem2 index \
   -p $INDEXDIR/GRCh38 \
   $GENOME