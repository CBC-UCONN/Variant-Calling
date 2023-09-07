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
module load singularity/3.10.0
module load vcflib/1.0.0-rc1 
module load bcftools/1.12
module load htslib/1.12

# input/output files, directories
INDIR=../../results/05_VariantCalling_longread/variants_clair3/reformatted

OUTDIR=../../results/05_VariantCalling_longread/variants_clair3_glnexus
mkdir -p $OUTDIR

# tell singularity where to write tmp files
export TMPDIR=/scratch/$USER
export SINGULARITY_TMPDIR=/scratch/$USER
mkdir -p $TMPDIR

# Pull the docker image. 
# singularity pull --dir $OUTDIR docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1
GLNEXUS_SIF=/isg/shared/databases/nfx_singularity_cache/glnexus_v1.4.1.sif

# GLNexus requires a configuration file, obtained from here: http://www.bio8.cs.hku.hk/clair3_trio/config/clair3.yml
CONFIG=clair3.yml


# run glnexus
singularity exec --bind /usr/lib/locale/ --bind $(pwd) $GLNEXUS_SIF glnexus_cli \
--dir $OUTDIR/scratch \
--config $CONFIG \
$INDIR/son.g.vcf.gz \
$INDIR/mom.g.vcf.gz \
$INDIR/dad.g.vcf.gz >$OUTDIR/joint.bcf

bcftools view -Oz $OUTDIR/joint.bcf >$OUTDIR/ashtrio_clair3.vcf.gz
