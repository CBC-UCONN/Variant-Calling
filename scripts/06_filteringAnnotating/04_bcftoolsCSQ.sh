#!/bin/bash
#SBATCH --job-name=bcftoolsCSQ
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load htslib/1.7
module load bcftools/1.20

# make a directory if it doesn't exist
INDIR=../../results/05_variantCalling/freebayes/
OUTDIR=../../results/06_Annotate/bcftoolsCSQ
mkdir -p ${OUTDIR}

GENOME=../../genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
VCFIN=../../results/05_variantCalling/freebayes/freebayes_filtered.vcf.gz
VCFOUT=${OUTDIR}/freebayes_annotated.vcf.gz

# first we need to download and fix up the GFF3 file. bcftools csq only works with ensembl GFF3 files (per the docs)
    # https://samtools.github.io/bcftools/howtos/csq-calling.html
# we'll grab it from here: https://useast.ensembl.org/Homo_sapiens/Info/Index
# and because it's an option, we'll only get chr20

wget -P ${OUTDIR} https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.chromosome.20.gff3.gz
gunzip ${OUTDIR}/Homo_sapiens.GRCh38.113.chromosome.20.gff3.gz
# fix up chromosome 20 names
sed -i 's/^20/chr20/' ${OUTDIR}/Homo_sapiens.GRCh38.113.chromosome.20.gff3

GFF=${OUTDIR}/Homo_sapiens.GRCh38.113.chromosome.20.gff3

# run bcftools csq
bcftools csq --phase a -f ${GENOME} -g ${GFF} ${VCFIN} -Oz -o ${VCFOUT}