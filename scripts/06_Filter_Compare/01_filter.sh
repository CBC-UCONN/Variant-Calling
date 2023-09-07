#!/bin/bash
#SBATCH --job-name=filter_vcfs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools/1.12
module load htslib/1.12
module load bedtools/2.29.0
module load GATK/4.3.0.0

INDIR=../../results/04_VariantCalling_shortread/
INDIRLONG=../../results/05_VariantCalling_longread/

OUTDIR=../../results/06_Filter_Compare/filtered_vcfs
mkdir -p $OUTDIR

TARGETS=../../results/03_AlignmentAndCoverage/short_read_coverage/targets.bed.gz

# using a subset (NC_000020.11:31400000-34400000) which excludes the most problematic areas

# jointly called variants for illumina data -------------------------

bcftools view -r NC_000020.11:31400000-34400000 $INDIR/variants_freebayes/ashtrio_fb.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/fb_ill_filter.vcf.gz

bcftools view -r NC_000020.11:31400000-34400000 $INDIR/variants_gatk/ashtrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/gatk_ill_filter.vcf.gz

# additionally filter bcftools output to match targets used for gatk and freebayes
bcftools view -r NC_000020.11:31400000-34400000 $INDIR/variants_bcftools/ashtrio.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bedtools intersect -header -wa -a stdin -b <(zcat $TARGETS) | bgzip -c > $OUTDIR/bcf_ill_filter.vcf.gz

bcftools view -r NC_000020.11:31400000-34400000 $INDIRLONG/variants_clair3_glnexus/ashtrio_clair3.vcf.gz  | bcftools filter -s LowQual -e '%QUAL<15' | bgzip -c > $OUTDIR/clair3_ONT_filter.vcf.gz

# glnexus apparently gets rid of the sequence dictionary. we have to add it back in. 
gatk UpdateVCFSequenceDictionary \
     -V $OUTDIR/clair3_ONT_filter.vcf.gz \
     --source-dictionary ../../results/03_AlignmentAndCoverage/bwa_align/son.bam \
     --output $OUTDIR/clair3_ONT_filter_dict.vcf.gz \
     --replace true

# son only to compare long reads against short reads for a single sample -------------------------

bcftools view -r NC_000020.11:31400000-34400000 $INDIRLONG/variants_clair3/son/merge_output.vcf.gz | bedtools intersect -header -wa -a stdin -b <(zcat $TARGETS) | bgzip -c > $OUTDIR/son_clair3_ont_filter.vcf.gz

# clair3's weird handling of case becomes a problem, so fix that...
zcat $OUTDIR/son_clair3_ont_filter.vcf.gz | \
	awk -F'\t' -v OFS='\t' '/^[^#]/{$4=toupper($4)} {$5=toupper($5)} {print $0}' | \
	bgzip >$OUTDIR/son_clair3_ont_filter.vcf.gz_casefixed && mv $OUTDIR/son_clair3_ont_filter.vcf.gz_casefixed $OUTDIR/son_clair3_ont_filter.vcf.gz

bcftools view -r NC_000020.11:31400000-34400000 $INDIR/variants_gatk_singlesample/son.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bedtools intersect -header -wa -a stdin -b <(zcat $TARGETS) | bgzip -c > $OUTDIR/son_gatk_ill_filter.vcf.gz

# index vcfs -------------------

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done
