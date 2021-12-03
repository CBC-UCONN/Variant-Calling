#!/bin/bash
#SBATCH --job-name=snpEff
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
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
OUTDIR=../annotated_vcfs
mkdir -p $OUTDIR

cd $OUTDIR

### functional prediction annotations with SnpEff

VCF=../filtered_vcfs/fb_vap.vcf.gz

# here -dataDir creates a directory where the hg38 database will be downloaded to
# the default directory cannot be written to by users

java -Xmx8G -jar $SNPEFF eff -dataDir $(pwd)/$OUTDIR/snpeff_data hg38 $VCF | bgzip -c >$OUTDIR/fb_vap.ann.vcf.gz
	
tabix -p vcf $OUTDIR/fb_vap.ann.vcf.gz
