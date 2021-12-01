#!/bin/bash 
#SBATCH --job-name=gatk_DBimport
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



hostname
date

# load required software
module load GATK/4.0

INDIR=../results/variants_gatk
OUTDIR=../results/variants_genomicsdb

#IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!
gatk --java-options "-Xmx10g -Xms4g" GenomicsDBImport \
  -V $INDIR/mom.g.vcf \
  -V $INDIR/dad.g.vcf \
  -V $INDIR/son.g.vcf \
  --genomicsdb-workspace-path $OUTDIR \
  --overwrite-existing-genomicsdb-workspace true \
  -L chr20:29400000-34400000

date