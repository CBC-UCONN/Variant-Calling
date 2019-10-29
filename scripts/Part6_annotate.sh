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

module load bcftools/1.6
module load htslib
module load snpEff/4.3q

# make a directory if it doesn't exist
mkdir -p annotated_vcfs
cd annotated_vcfs


### functional prediction annotations with SnpEff

VCF=../filtered_vcfs/fb_filter.vcf.gz

# here -dataDir creates a directory where the hg38 database will be downloaded to
# the default directory cannot be written to by users

java -Xmx8G -jar $SNPEFF eff -dataDir ~/cbc_projects/vc_workshop/annotated_vcfs/temp_data hg38 $VCF > fb_filter.ann.vcf.gz
	





### Annotating variants with dbSNP rsids. 


# get the dbsnp set for chromosome 20
	# do a bunch of reformatting to make sure it works with the following steps
	
	# download only a section of chr20 from dbsnp
	tabix -h ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz 20:28000000-35000000 | \
	sed 's/^20/chr20/' | \
	bgzip -c >chr20.dbsnp.vcf.gz
	# update the sequence dictionary
	gatk UpdateVCFSequenceDictionary -V chr20.dbsnp.vcf.gz --source-dictionary ../align_pipe/son.bam --output chr20.dbsnp.contig.vcf.gz --replace=true
	# make sure all the filter columns say PASS so we can filter below and not lose all the variants
	bcftools filter -s LowQual chr20.dbsnp.contig.vcf.gz | bgzip -c >chr20.dbsnp.filter.vcf.gz
	# index
	tabix -p vcf chr20.dbsnp.filter.vcf.gz
	# remove intermediate files
	rm chr20.dbsnp.contig.vcf.gz chr20.dbsnp.vcf.gz
