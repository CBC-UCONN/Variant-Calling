#!/bin/bash 
#SBATCH --job-name=clair3_gvcf
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load htslib/1.12
module load singularity/3.10.0

# we're going to run clair3 inside of a software container. 
# we already have the container on the xanadu cluster in a shared location
# to get your own copy of the container, load the singularity module and then:
 # singularity pull docker://hkubal/clair3:latest
 # that will result in a container called "clair3_latest.sif" which we use below


# Set the number of CPUs to use
THREADS="4"

# input/output files, directories
INDIR=../../results/03_AlignmentAndCoverage/minimap2/

OUTDIR=../../results/05_VariantCalling_longread/variants_clair3
mkdir -p $OUTDIR

GENOME=../../genome/GRCh38_latest_genomic.fna

CLAIR3_CONTAINER=/isg/shared/databases/nfx_singularity_cache/clair3_latest.sif

TARGETS=$OUTDIR/targets.bed
  zcat ../../results/03_AlignmentAndCoverage/short_read_coverage/targets.bed.gz >$TARGETS

# run clair3

# son
SAM=son
singularity exec $CLAIR3_CONTAINER \
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=/opt/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM \
  --bed_fn=$TARGETS

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

# mom
SAM=mom
singularity exec $CLAIR3_CONTAINER \
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=/opt/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM \
  --bed_fn=$TARGETS

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

# dad
SAM=dad
singularity exec $CLAIR3_CONTAINER \
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=/opt/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM \
  --bed_fn=$TARGETS

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

date

