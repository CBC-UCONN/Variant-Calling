#!/bin/bash 
#SBATCH --job-name=coverage
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



hostname
date

# load required software

module load bedtools
module load bamtools
module load htslib
module load samtools

# define and/or create input, output directories

INDIR=../align_pipe
OUTDIR=../coverage_stats
mkdir -p $OUTDIR

# genome index file from samtools faidx
FAI=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta.fai

# make a "genome" file, required by bedtools makewindows command, set variable for location
GFILE=$OUTDIR/Homo_sapiens_assembly38.fasta.genome
cut -f 1-2 $FAI > $GFILE

# make 1kb window bed file, set variable for location
WIN1KB=$OUTDIR/hg38_1kb.bed
bedtools makewindows -g $GFILE -w 1000 >$WIN1KB


# make a list of bam files
find $INDIR -name "*bam" >$OUTDIR/bam.list

# summarize coverage as the number of fragments mapping to 1kb windows across the genome
# pipe:
	# 1) merge bam files
	# 2) filter by quality and proper pairing
	# 3) convert alignments to bed format
	# 4) map alignments to 1kb windows, counting (but also getting the mean and median of the mapping quality score from column 5)

bamtools merge -list $OUTDIR/bam.list | \
bamtools filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
bedtools bamtobed -i stdin | \
bedtools map \
-a $WIN1KB \
-b stdin \
-c 5 -o mean,median,count \
-g $GFILE | \
bgzip >$OUTDIR/coverage_1kb.bed.gz

# bgzip compress and tabix index the resulting file
tabix -p bed $OUTDIR/coverage_1kb.bed.gz

# select and merge outlier windows (after deciding what is an outlier by looking at the distribution in R)
zcat $OUTDIR/coverage_1kb.bed.gz | awk '$6 < 850 || $6 > 2550' | bedtools merge | bgzip >$OUTDIR/coverage_outliers.bed.gz 
tabix -p bed $OUTDIR/coverage_outliers.bed.gz

# select and merge target windows (inverse of "outlier windows" above)
zcat $OUTDIR/coverage_1kb.bed.gz | awk '$6 > 850 && $6 < 2550' | bedtools merge | bgzip >$OUTDIR/targets.bed.gz 
tabix -p bed $OUTDIR/targets.bed.gz


# calculate per-base coverage as well	
bamtools merge -list $OUTDIR/bam.list | \
bamtools filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | 
samtools depth -d 20000 /dev/stdin | \
bgzip >$OUTDIR/per_base_coverage.txt.gz

date




