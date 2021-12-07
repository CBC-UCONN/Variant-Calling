#!/bin/bash 
#SBATCH --job-name=glnexus_clair3
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
module load singularity/3.7.1
module load vcflib/1.0.0-rc1 
module load bcftools/1.12
module load htslib/1.12

# input/output files, directories
INDIR=../results/variants_clair3_gvcf

OUTDIR=../results/variants_clair3_joint
mkdir -p $OUTDIR

# tell singularity where to write tmp files
export TMPDIR=/scratch/$USER
export SINGULARITY_TMPDIR=/scratch/$USER
mkdir -p $TMPDIR

# Pull the docker image
singularity pull --dir $OUTDIR docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1

# run glnexus
singularity exec --bind /usr/lib/locale/ --bind $(pwd) $OUTDIR/glnexus_v1.4.1.sif glnexus_cli \
--dir $OUTDIR/scratch \
--config gatk \
$INDIR/son.g.vcf.gz \
$INDIR/mom.g.vcf.gz \
$INDIR/dad.g.vcf.gz >$OUTDIR/joint.bcf

bcftools view -Oz $OUTDIR/joint.bcf >$OUTDIR/ashtrio_clair3.vcf.gz