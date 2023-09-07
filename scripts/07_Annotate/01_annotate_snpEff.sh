#!/bin/bash
#SBATCH --job-name=snpEff
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
module load bcftools/1.6
module load snpEff/4.3q

# make a directory if it doesn't exist
OUTDIR=../../results/07_Annotate/
mkdir -p $OUTDIR

cd $OUTDIR

### functional prediction annotations with SnpEff

VCFIN=../06_Filter_Compare/filtered_vcfs/fb_vap.vcf.gz
VCF=fb_vap_rename.vcf.gz
VCFANN=fb_vap.ann.vcf.gz

# we're using the refseq version of the genome, and the chromosomes are all referred to by refseq accession numbers, e.g. NC_00001.11, instead of chr1
# so we need to fix that. we're only working with chromosome 20 in this tutorial, so it's easy. 

zcat $VCFIN | sed 's/NC_000020.11/chr20/' | bgzip >$VCF

# here -dataDir creates a directory where the hg38 database will be downloaded to
# the default directory cannot be written to by users

java -Xmx8G -jar $SNPEFF eff -dataDir $(pwd)/$OUTDIR/snpeff_data hg38 $VCF | bgzip -c >$VCFANN
	
tabix -p vcf $VCFANN
