#!/bin/bash
#SBATCH --job-name=summarize_vcfs
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
module load bcftools/1.12
module load htslib/1.12
module load vcflib/1.0.0-rc1

cd ../../results/06_Filter_Compare/filtered_vcfs

##############################
# summarize variant call sets:
##############################

# use VT peek to tally up variant types
vt peek fb_ill_filter.vcf.gz
vt peek -f PASS fb_ill_filter.vcf.gz
vt peek -f PASS gatk_ill_filter.vcf.gz
vt peek -f PASS bcf_ill_filter.vcf.gz
vt peek -f PASS clair3_ONT_filter.vcf.gz

vt peek -f PASS son_clair3_ont_filter.vcf.gz  
vt peek -f PASS son_gatk_ill_filter.vcf.gz

# use VT profile_mendelian to look at mendelian error rates in the trio
# first create a pedigree file:
echo -e 'giabct\tson\tdad\tmom\t1\t-9' >ct.ped
echo -e 'giabct\tmom\t0\t0\t2\t-9' >>ct.ped
echo -e 'giabct\tdad\t0\t0\t1\t-9' >>ct.ped

vt profile_mendelian -f PASS -p ct.ped fb_ill_filter.vcf.gz
vt profile_mendelian -f PASS -p ct.ped gatk_ill_filter.vcf.gz
vt profile_mendelian -f PASS -p ct.ped bcf_ill_filter.vcf.gz
vt profile_mendelian -f PASS -p ct.ped clair3_ONT_filter.vcf.gz

