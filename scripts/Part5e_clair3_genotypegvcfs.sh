#!/bin/bash 
#SBATCH --job-name=clair3_gatk_genotype
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
module load htslib/1.12

# input/output directories
INDIR=../results/variants_clair3_gvcf

DBDIR=../results/variants_clair3_genomicsdb
mkdir -p $DBDIR

OUTDIR=../results/variants_clair3_joint_gatk
mkdir -p $OUTDIR

# import databases
#IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!
gatk --java-options "-Xmx10g -Xms4g" GenomicsDBImport \
  -V $INDIR/mom.g.vcf.gz \
  -V $INDIR/dad.g.vcf.gz \
  -V $INDIR/son.g.vcf.gz \
  --genomicsdb-workspace-path $DBDIR \
  --overwrite-existing-genomicsdb-workspace true \
  -L chr20:29400000-34400000


# do genotyping
GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

gatk GenotypeGVCFs \
    -R $GEN \
    -V gendb://$DBDIR \
    -O $OUTDIR/ashtrio.vcf 

bgzip $OUTDIR/ashtrio.vcf 
tabix -p vcf $OUTDIR/ashtrio.vcf.gz

date