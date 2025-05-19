#!/bin/bash 
#SBATCH --job-name=download_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date


OUTDIR=../../genome
mkdir -p $OUTDIR
cd $OUTDIR

# this is modified version of GRCh38 that corrects some errors. GIAB uses this for their benchmarking work. 
    # "with masked false duplications and contaminations, as well as decoy sequences from CHM13, which we are now using for GIAB analyses"
    #  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02863-7
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
gunzip *gz