# some random bash lines you can use to dig into VCF files
    # on large VCF files some lines should be done as batch jobs

# load bcftools
module load bcftools/1.20

# make some variables to point at our VCF files
BCFTOOLS=../../results/05_variantCalling/bcftools/bcftools.vcf.gz
FREEBAYES=../../results/05_variantCalling/freebayes/freebayes.vcf.gz
GATK=../../results/05_variantCalling/gatk/gatk.vcf.gz

######## look at the contents of a VCF file with `bcftools view` ###########

# print the header and the records
bcftools view $BCFTOOLS

# print only the header with -h
bcftools view -h $BCFTOOLS

# suppress the header with -H and print only the records
bcftools view -H $BCFTOOLS | head

# see records from a region:
bcftools view -H $BCFTOOLS chr20:31000000-31010000

############## summarize the vcf with `bcftools stats` ##############
    # when you have a large VCF, run this in a batch script
    
# pull out the summary numbers
bcftools stats $BCFTOOLS | grep "^SN"

# pull out the ts/tv ratio
bcftools stats $BCFTOOLS | grep "TSTV"

# filter out QUAL < 30
    # filter flags can be applied to most bcftools commands
bcftools stats -e "QUAL < 30" $BCFTOOLS | grep "TSTV"

############### extract particular bits of information into a table with `bcftools query` #############

bcftools query \
    --print-header \
    -r chr20:33000000- \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%GT\t]' \
    $FREEBAYES | 
    head


bcftools query \
    --print-header \
    -r chr20:33000000- \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AB\t[%GT\t][%DP\t]' \
    $FREEBAYES | \
    head

############## calculating the frequency of mendelian violations #####################

bcftools +mendelian2 -r chr20:29400000-34400000 -p son,dad,mom -m c $FREEBAYES