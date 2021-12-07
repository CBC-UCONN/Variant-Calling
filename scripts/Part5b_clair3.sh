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
source ~/.bashrc
conda activate clair3

# Set the number of CPUs to use
THREADS="4"

# input/output files, directories
INDIR=../results/align_ont

OUTDIR=../results/variants_clair3_gvcf
mkdir -p $OUTDIR

GENOME=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

# run clair3

# son
SAM=son
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=$CONDA_PREFIX/bin/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

# mom
SAM=mom
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=$CONDA_PREFIX/bin/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

# dad
SAM=dad
run_clair3.sh \
  --bam_fn=$INDIR/$SAM.bam \
  --ref_fn=$GENOME \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=$CONDA_PREFIX/bin/models/r941_prom_hac_g360+g422 \
  --output=$OUTDIR/$SAM \
  --gvcf \
  --sample_name=$SAM

cp $OUTDIR/$SAM/merge_output.gvcf.gz $OUTDIR/$SAM.g.vcf.gz
tabix -p vcf $OUTDIR/$SAM.g.vcf.gz

date











