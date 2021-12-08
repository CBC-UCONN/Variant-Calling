#!/bin/bash 
#SBATCH --job-name=pepper
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load singularity/3.7.1

# Set the number of CPUs to use
THREADS="1"

# input/output files, directories
INDIR=../results/align_ont

OUTDIR=../results/variants_pepper
mkdir -p $OUTDIR

GENOME=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

# set up singularity
export TMPDIR=/scratch/$USER
export SINGULARITY_TMPDIR=/scratch/$USER
mkdir -p $TMPDIR

# Pull the docker images
singularity pull --dir $OUTDIR docker://kishwars/pepper_deepvariant:r0.5

# PEPPER-Margin-DeepVariant

# son
SAM=son
singularity exec --bind /usr/lib/locale/ --bind $(pwd) $OUTDIR/pepper_deepvariant_r0.5.sif \
run_pepper_margin_deepvariant call_variant \
-b $INDIR/$SAM.bam \
-f $GENOME \
-t $THREADS \
--ont \
-o $OUTDIR/$SAM \
-r chr20:29400000-34400000 \
-p $SAM \
-s $SAM

# mom
SAM=mom
singularity exec --bind /usr/lib/locale/ --bind $(pwd) $OUTDIR/pepper_deepvariant_r0.5.sif \
run_pepper_margin_deepvariant call_variant \
-b $INDIR/$SAM.bam \
-f $GENOME \
-t $THREADS \
--ont \
-o $OUTDIR/$SAM \
-r chr20:29400000-34400000 \
-p $SAM \
-s $SAM

# dad
SAM=dad
singularity exec --bind /usr/lib/locale/ --bind $(pwd) $OUTDIR/pepper_deepvariant_r0.5.sif \
run_pepper_margin_deepvariant call_variant \
-b $INDIR/$SAM.bam \
-f $GENOME \
-t $THREADS \
--ont \
-o $OUTDIR/$SAM \
-r chr20:29400000-34400000 \
-p $SAM \
-s $SAM
