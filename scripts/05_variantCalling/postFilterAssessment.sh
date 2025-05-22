# load bcftools 
module load bcftools/1.20

# set variables
BCFTOOLS=../../results/05_variantCalling/bcftools/bcftools.vcf.gz
FREEBAYES=../../results/05_variantCalling/freebayes/freebayes.vcf.gz
GATK=../../results/05_variantCalling/gatk/gatk.vcf.gz

BCFILTER=../../results/05_variantCalling/bcftools/bcftools_filtered.vcf.gz
FBFILTER=../../results/05_variantCalling/freebayes/freebayes_filtered.vcf.gz
GATKFILTER=../../results/05_variantCalling/gatk/gatk_filtered.vcf.gz

# ts/tv
# unfiltered
bcftools stats $FREEBAYES | grep "TSTV"
# filtered
bcftools stats $FBFILTER | grep "TSTV"

# mendelian violations
# unfiltered
bcftools +mendelian2 -p son,dad,mom -m c $FREEBAYES
# filtered
bcftools +mendelian2 -p son,dad,mom -m c $FBFILTER

# genetic diversity, theta. expected value is ~ 0.001.
    # use Watterson estimator: https://en.wikipedia.org/wiki/Watterson_estimator
# this is a very, very rough approach
    # theta = (S / alpha) / nsites
    # alpha = sum(1/(1:n-1))
        # S is number of variable sites
        # n is number of sampled alleles
            # here we have 3 individuals, but only 2 independent ones, so say 4 alleles. 
        # nsites is the number of sites assayed for variation. 

# unfiltered
    # nsites=5000001
    # alpha = 1.833
# for S:
bcftools stats $FREEBAYES | grep "^SN" | grep records
    # S = 71579
# theta = (71579 / 1.833) / 5000001 = 0.007810037

# filtered
# alpha = 1.833
# nsites:
TARGETS=../../results/04_alignQC/coverage/targets.bed.gz
zcat $TARGETS | awk '{n+=$3-$2}END{print n}'
    # nsites = 4416672
# for S:
bcftools stats $FBFILTER | grep "^SN" | grep records
    # S = 7718
# theta = (7718 / 1.833) / 4416672 = 0.0009533386