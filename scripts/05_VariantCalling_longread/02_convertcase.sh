#!/bin/bash 
#SBATCH --job-name=clair3_reformat
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load singularity/3.10.0
module load vcflib/1.0.0-rc1 
module load bcftools/1.12
module load htslib/1.12


# directories, files
INDIR=../../results/05_VariantCalling_longread/variants_clair3

OUTDIR=../../results/05_VariantCalling_longread/variants_clair3/reformatted
mkdir -p $OUTDIR

# change lower case letters to upper case
SAM=son
zcat $INDIR/$SAM.g.vcf.gz | \
	awk -F'\t' -v OFS='\t' '/^[^#]/{$4=toupper($4)} {$5=toupper($5)} {print $0}' | \
	bgzip >${OUTDIR}/$SAM.g.vcf.gz

SAM=dad
zcat $INDIR/$SAM.g.vcf.gz | \
	awk -F'\t' -v OFS='\t' '/^[^#]/{$4=toupper($4)} {$5=toupper($5)} {print $0}' | \
	bgzip >${OUTDIR}/$SAM.g.vcf.gz

SAM=mom
zcat $INDIR/$SAM.g.vcf.gz | \
	awk -F'\t' -v OFS='\t' '/^[^#]/{$4=toupper($4)} {$5=toupper($5)} {print $0}' | \
	bgzip >${OUTDIR}/$SAM.g.vcf.gz

for file in $(ls $OUTDIR/*vcf.gz); do tabix -p vcf $file; done
