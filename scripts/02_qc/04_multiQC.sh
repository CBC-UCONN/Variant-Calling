#!/bin/bash
#SBATCH --job-name=multiqc_raw
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#######################
# run multiqc
#######################

# aggregate fastqc reports

# input, output directories
INDIR=../../results/02_qc/

OUTDIR=../../results/02_qc/multiqc
mkdir -p $OUTDIR


module load MultiQC/1.9

multiqc --outdir $OUTDIR $INDIR