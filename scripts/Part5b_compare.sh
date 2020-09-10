#!/bin/bash
#SBATCH --job-name=bcf_vcf
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load vt/0.57721
module load bcftools/1.6
module load htslib
module load vcflib

##############################
# check out filtering results:
##############################

# use VT peek
vt peek fb_filter.vcf.gz
vt peek -f PASS fb_filter.vcf.gz
vt peek -f PASS gatk_filter.vcf.gz
vt peek -f PASS bcf_filter.vcf.gz

# use VT profile_mendelian
# first create a pedigree file:
echo giabct	son	dad	mom	1 -9 >ct.ped
echo giabct	son	0	0	1 -9 >>ct.ped
echo giabct	son	0	0	1 -9 >>ct.ped

vt profile_mendelian -f PASS -p ct.ped fb_filter.vcf.gz
vt profile_mendelian -f PASS -p ct.ped gatk_filter.vcf.gz
vt profile_mendelian -f PASS -p ct.ped bcf_filter.vcf.gz

##############################
# compare variant sets
##############################
# for two sets: VT partition

vt partition -f PASS fb_filter.vcf.gz gatk_filter.vcf.gz
vt partition -f PASS fb_filter.vcf.gz bcf_filter.vcf.gz
vt partition -f PASS gatk_filter.vcf.gz bcf_filter.vcf.gz

# for three sets: VT multipartition

vt multi_partition -f PASS fb_filter.vcf.gz gatk_filter.vcf.gz bcf_filter.vcf.gz

# concordance is low because many complex variants and small haplotypes are represented differently
	#  bcftools norm and vcflib's vcfallelicprimitives can help standardize the representation

GEN=/UCHC/PublicShare/CBC_Tutorials/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

bcftools norm -f $GEN fb_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >fb_vap.vcf.gz
bcftools norm -f $GEN gatk_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >gatk_vap.vcf.gz
bcftools norm -f $GEN bcf_filter.vcf.gz | vcfallelicprimitives | vcfstreamsort | bgzip >bcf_vap.vcf.gz

for file in *vap.vcf.gz
do tabix -f -p vcf $file
done

# repeat vt partition and multipartition
vt partition -f PASS fb_vap.vcf.gz gatk_vap.vcf.gz
vt partition -f PASS fb_vap.vcf.gz bcf_vap.vcf.gz
vt partition -f PASS gatk_vap.vcf.gz bcf_vap.vcf.gz

vt multi_partition -f PASS fb_vap.vcf.gz gatk_vap.vcf.gz bcf_vap.vcf.gz


# extract intersections with bcftools isec
bcftools isec -f PASS -p isec_fb_gatk fb_vap.vcf.gz gatk_vap.vcf.gz

bcftools view -H isec_fb_gatk/0000.vcf | head -n 10

# look at a few incongruent markers from isec_fb_gatk/0000.vcf, which contains variants called by fb but not gatk

# chr20	31577045	.	T	C
	# called by fb, bcf, but not gatk. no clue why. 
	# apparently lots of support for C alternate allele. 
	# heterozygous parent has RO:35, AO:21





