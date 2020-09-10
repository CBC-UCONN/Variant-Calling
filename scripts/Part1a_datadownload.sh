#!/bin/bash 
#SBATCH --job-name=download_data
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load samtools
module load bedtools

# script assumes it is run in directory vc_workshop/scripts

OUTDIR=../rawdata/
mkdir -p $OUTDIR

# download a subregion of chinese GIAB trio: chr20:29400000-34400000
# sort reads by name, convert to fastq. 
# there will be lots of errors from bedtools resulting from discordant reads. 
	# this is not an issue in this specific tutorial case, but could be worrisome in other contexts. 
# note: the son is aligned to a different reference genome, so the coordinates don't quite line up. 

# son
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'
samtools view -uh $SON chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/son.1.fq -fq2 $OUTDIR/son.2.fq

# mom
MOM='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam'
samtools view -uh $MOM chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/mom.1.fq -fq2 $OUTDIR/mom.2.fq

# dad
DAD='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam'
samtools view -uh $DAD chr20:29400000-34400000 | \
samtools sort -n - | \
bedtools bamtofastq -i /dev/stdin/ -fq $OUTDIR/dad.1.fq -fq2 $OUTDIR/dad.2.fq

# get rid of bam indexes that were also downloaded
rm $OUTDIR/*bai
rm *bam.bai
