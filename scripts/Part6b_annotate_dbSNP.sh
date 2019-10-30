#!/bin/bash
#SBATCH --job-name=snpEff
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load htslib/1.7
module load bcftools/1.6
module load GATK/4.1.3.0

### Annotating variants with dbSNP rsids. 

# get the dbsnp set for chromosome 20
	
	# download only a section of chr20 from dbsnp
	tabix -h ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz 20:28000000-35000000 | \
	sed 's/^20/chr20/' | \
	bgzip -c >chr20.dbsnp.vcf.gz
	tabix -p vcf chr20.dbsnp.vcf.gz
	# update the sequence dictionary
	gatk UpdateVCFSequenceDictionary -V chr20.dbsnp.vcf.gz --source-dictionary ../align_pipe/son.bam --output chr20.dbsnp.contig.vcf.gz --replace=true
	tabix -p vcf chr20.dbsnp.contig.vcf.gz
	# remove intermediate files
	rm chr20.dbsnp.vcf.gz*

# annotate:
tabix -p vcf fb_filter.ann.vcf.gz
bcftools annotate -c ID \
--output-type z \
-a chr20.dbsnp.contig.vcf.gz \
fb_filter.ann.vcf.gz >fb_filter.ann.RSID.vcf.gz


