#!/bin/bash
#SBATCH --job-name=compare_vcfs
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

module load vt/0.57721
module load bcftools/1.12
module load htslib/1.12
module load vcflib/1.0.0-rc1

cd ../../results/06_Filter_Compare/filtered_vcfs

##############################
# compare variant sets
##############################
# for two sets: VT partition

vt partition -f PASS fb_ill_filter.vcf.gz gatk_ill_filter.vcf.gz
vt partition -f PASS fb_ill_filter.vcf.gz bcf_ill_filter.vcf.gz
vt partition -f PASS fb_ill_filter.vcf.gz clair3_ONT_filter_dict.vcf.gz
vt partition -f PASS gatk_ill_filter.vcf.gz bcf_ill_filter.vcf.gz

vt partition -f PASS son_clair3_ont_filter.vcf.gz son_gatk_ill_filter.vcf.gz


# for three sets: VT multipartition

vt multi_partition -f PASS fb_ill_filter.vcf.gz gatk_ill_filter.vcf.gz bcf_ill_filter.vcf.gz

# concordance is low because many complex variants and small haplotypes are represented differently in the call sets. 
	#  bcftools norm and vcflib's vcfallelicprimitives can help standardize the representation

GEN=../../../genome/GRCh38_latest_genomic.fna

bcftools norm -f $GEN fb_ill_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >fb_vap.vcf.gz
bcftools norm -f $GEN gatk_ill_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >gatk_vap.vcf.gz
bcftools norm -f $GEN bcf_ill_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >bcf_vap.vcf.gz
bcftools norm -f $GEN clair3_ONT_filter_dict.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >clair3_ONT_vap.vcf.gz
bcftools norm -f $GEN son_clair3_ont_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >son_clair3_vap.vcf.gz
bcftools norm -f $GEN son_gatk_ill_filter.vcf.gz  | vcfallelicprimitives | vcfstreamsort | bgzip >son_gatk_vap.vcf.gz


for file in *vap.vcf.gz
do tabix -f -p vcf $file
done

# repeat vt partition and multipartition
vt partition -f PASS fb_vap.vcf.gz gatk_vap.vcf.gz
vt partition -f PASS fb_vap.vcf.gz bcf_vap.vcf.gz
vt partition -f PASS fb_vap.vcf.gz clair3_ONT_vap.vcf.gz
vt partition -f PASS gatk_vap.vcf.gz bcf_vap.vcf.gz

vt multi_partition -f PASS fb_vap.vcf.gz gatk_vap.vcf.gz bcf_vap.vcf.gz

vt multi_partition -f PASS fb_vap.vcf.gz gatk_vap.vcf.gz clair3_ONT_vap.vcf.gz

# extract intersections with bcftools isec
bcftools isec -f PASS -p isec_fb_clair3 fb_vap.vcf.gz clair3_ONT_vap.vcf.gz

bcftools view -H isec_fb_clair3/0000.vcf | head -n 10 | cut -f 1-6
bcftools view -H isec_fb_clair3/0001.vcf | head -n 10 | cut -f 1-6
# files produced by isec:

# isec_fb_clair3/0000.vcf   for records private to  fb_vap.vcf.gz
# isec_fb_clair3/0001.vcf   for records private to  clair3_ONT_vap.vcf.gz
# isec_fb_clair3/0002.vcf   for records from fb_vap.vcf.gz shared by both   fb_vap.vcf.gz clair3_ONT_vap.vcf.gz
# isec_fb_clair3/0003.vcf   for records from clair3_ONT_vap.vcf.gz shared by both fb_vap.vcf.gz clair3_ONT_vap.vcf.gz

# look at a few incongruent markers from isec_fb_clair3/0000.vcf, which contains variants called by fb but not gatk


# now let's look at some of the long-read data and single-sample (gatk) short read variants

vt peek -f PASS clair3_son_ont_filter.vcf.gz
vt peek -f PASS pepper_son_ont_filter.vcf.gz
vt peek -f PASS gatk_son_ont_filter.vcf.gz

vt partition -f "PASS && VTYPE==SNP" clair3_son_ont_filter.vcf.gz pepper_son_ont_filter.vcf.gz
vt partition -f "PASS && VTYPE==SNP" clair3_son_ont_filter.vcf.gz gatk_son_ont_filter.vcf.gz
vt partition -f "PASS && VTYPE==SNP" pepper_son_ont_filter.vcf.gz gatk_son_ont_filter.vcf.gz

vt multi_partition -f "PASS && VTYPE==SNP" clair3_son_ont_filter.vcf.gz pepper_son_ont_filter.vcf.gz gatk_son_ont_filter.vcf.gz

# extract intersections with bcftools isec
bcftools isec -i "FILTER='PASS' && TYPE='SNP'" -p isec_clair3_pepper clair3_son_ont_filter.vcf.gz pepper_son_ont_filter.vcf.gz

bcftools isec -i "FILTER='PASS' && TYPE='SNP'" -p isec_clair3_gatk clair3_son_ont_filter.vcf.gz gatk_son_ont_filter.vcf.gz

# have a look at some variants called by clair3/ont, but not gatk/illumina
bcftools view -H isec_clair3_gatk/0000.vcf | head -n 10
