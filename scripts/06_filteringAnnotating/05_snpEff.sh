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

module load htslib/1.16
module load bcftools/1.20
module load snpEff/4.3q

# make a directory if it doesn't exist
OUTDIR=../../results/06_Annotate/snpEff
mkdir -p $OUTDIR

# easier to run snpEff from the work directory
cd $OUTDIR

### functional prediction annotations with SnpEff

#vcf in/out
VCFIN=../../05_variantCalling/freebayes/freebayes_filtered.vcf.gz
VCFANN=freebayes_annotated.vcf.gz

# here -dataDir creates a directory where the hg38 database will be downloaded to
# the default directory cannot be written to by users
# apparently snpEff wants an absolute file path. 

java -Xmx8G -jar ${SNPEFF} eff -dataDir $(pwd)/snpeff_data hg38 $VCFIN | bgzip -c >${VCFANN}
	
tabix -p vcf ${VCFANN}