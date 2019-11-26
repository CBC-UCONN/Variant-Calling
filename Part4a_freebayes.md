# Variant discovery using freebayes

## Introduction

This section of the tutorial introduces variant calling using `freebayes`. It assumes you have completed Part 1 and Part 3 and have QC'ed, aligned, post-processed bam files in the directory structure first created in Part 1. 

Steps here will use the following software packages:

- [ freebayes ](https://github.com/ekg/freebayes)
- [ bedtools ](https://bedtools.readthedocs.io/en/latest/)
- [ bamtools ](https://github.com/pezmaster31/bamtools)
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)
- [ tabix ](http://www.htslib.org/doc/tabix.html)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Decide where to call variants ](#Decide-where-to-call-variants)
3.    [ Call variants ](#Call-variants)


## Motivation

In [Part 2](/Part2_bcftools.md) we used `bcftools` to call variants. You will remember the two step protocol in `bcftools` is to first summarize the read pileup, base by base, across the genome, and then evaluate the evidence for variation at each site. This means `bcftools` takes the alignments in the read pileup literally. 

Two factors combine to make this problematic for variant calling. First, indels and complex variants (indels + SNPs) may not have a single unambiguous representation. Second, we align reads (or read pairs) independently. This means that a single underlying biological sequence that differs from the reference may be represented by multiple different alignments in the read pileup. This can lead to both false negative and false positive variant calls. 

`bcftools` deals with this problem by downgrading base qualities to account for uncertainty in the sequence alignment (using ["base alignment qualities"](https://www.ncbi.nlm.nih.gov/pubmed/21320865)). This reduces the probability of false positive calls, but it also reduces the sensitivity to indels and nearby SNP variants. 

The variant caller we will use here, `freebayes` and another, the GATK's `HaplotypeCaller` (discussed in [Part 4b](Part4b_gatk.md)) take different approaches to this problem. Instead of calling variants at single sites, `freebayes` calls __haplotypes__. Haplotypes are combinations of alleles at different sites found on a single chromosome. To do this, `freebayes` takes in the aligned reads overlapping a short region of the genome, standardizes the representation of indels, identifies candidate haplotypes from the reads (which must be shorter than the read length), and then evaluates the evidence for variation at the site using a Bayesian model. 

INSERT FIGURE

The key innovation here is the use of haplotype alleles. These allow the evidence for indels and complex variants with many possible alignments to be evaluated more rigorously, improving the sensitivity and specificity for these calls. 

`freebayes` has a few other advantages as well. The biggest among them is that it is designed with piping in mind. It can accept a single stream of BAM input merged from multiple samples. That means that filtering or other alignment manipulations can be performed on the fly, delivering data to the caller without generating intermediate alignment files. We will see how this can work below. 


## Update your working directory

First we'll create a new directory to hold the results we'll generate. Make sure you're in the directory vc_workshop and type:

```bash
mkdir -p variants_freebayes coverage_stats
```

## Decide where to call variants

The most common approach to discovering genomic variation, short-read sequencing followed by reference mapping, presents a few challenges:
- Most genomes are drafts. This means that many regions present in the true genome of the reference individual are absent from the reference sequence. 
- Many genomes contain regions of repetitive or highly similar sequence. 
- The genomic composition of individuals within a species can be highly variable, with some regions copied, sometimes many times, or deleted. 
- Some regions present in the reference may simply be resistant to sequencing. 
- Finally, the sizes of fragments typically sequenced are often smaller than the size of repetitive sequence motifs. 

Because variant callers must assume (conditional on the mapping quality score) that reads aligned to the reference sequence were derived from the homologous region of another individual, these issues mean that your sequence alignments are likely to have problematic regions that lead to unreliable inference. 

These regions come in three main flavors:
- Regions that have more than expected sequencing coverage because they bear some sequence similarity to regions that are absent (or underrepresented) in the reference sequence. 
- Regions that have less than expected sequencing coverage because they have been deleted or are resistant to sequencing. 
- Regions that have the expected sequencing coverage, but reads mapping there have low mapping qualities because there are other very similar regions present in the reference sequence. 

In regions with excess coverage, true biological variation may be discovered, but it will be unclear if that variation represents polymorphism within an individual, divergence between copies within the same genome, or both. In any case, the true location of the discovered variation will also be unknown. Mapping quality scores may not be of any help here, as high sequence similarity to a missing region can result in high mapping quality scores. Furthermore, the excess pileup of reads in some of these regions can be extreme, with hundreds of times the expected coverage. This can dramatically slow the variant caller and balloon memory usage, unncessarily stalling an analysis. 

In regions with deficient coverage, sequence variation may be discovered, but called genotypes can be inaccurate. For example, a 50kb deletion will not be detected by a variant caller such as `freebayes`. If a diploid individual is heterozygous for that deletion, any SNP or short indel genotypes called within that deletion will be called as homozygous. 

Because we are likely to need to filter these regions out later, and the very high coverage ones may be computational bottlenecks, it is a good idea (though not always necessary depending on the application) to exclude them from initial variant calling. 

There are many ways we can accomplish this, here we'll just try to identify 1kb windows with aberrantly high and low coverage and exclude those. We'll use `bedtools`, a really excellent program for manipulating genomic windows, and bamtools, which can merge and filter bam files. 

```bash
# genome index file from samtools faidx
FAI=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta.fai

# make a "genome" file, per bedtools makewindows command
cut -f 1-2 $FAI > Homo_sapiens_assembly38.fasta.genome
GFILE=Homo_sapiens_assembly38.fasta.genome

# tell bedtools to output a "BED" file of 1kb non-overlapping windows
bedtools makewindows -g $GFILE -w 1000 >hg38_1kb.bed
```

An aside about BED format:

```bash
chr1	0	1000
chr1	1000	2000
chr1	2000	3000
```

The simplest BED formatted file is a tab delimited table giving genomic intervals. Column 1 gives the sequence, column 2 gives the start of the window (0-indexed) and column 3 gives the end of the window (1-indexed). The 0 vs 1 indexing switch between start and stop positions often trips people up. It means the first base in any given window is actually 'start+1', but the last base is simply 'end'. Following columns can contain anything. 

Now that we have a BED file with 3.2 million 1kb windows, (`wc -l hg38_1kb.bed`) we have to find out how much sequencing coverage they each have. 

```bash
# set variable for 1kb window folder
WIN1KB=../coverage_stats/hg38_1kb.bed

# go to sequence alignment folder
cd ../align_pipe

# make a list of bam files
ls *bam >bam.list

# pipe:
	# 1) merge bam files
	# 2) filter by quality and proper pairing
	# 3) convert alignments to bed format
	# 4) map alignments to 1kb windows, counting (but also getting the mean and median of the mapping quality score from column 5)

bamtools merge -list bam.list | \
bamtools filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
bedtools bamtobed -i stdin | \
bedtools map \
-a $WIN1KB \
-b stdin \
-c 5 -o mean,median,count \
-g $GFILE \
>../coverage_stats/coverage_1kb.bed
```

This pipeline merges the bam files into a single stream, passes them to a filter, which keeps only properly paired reads with mapping quality > 30, converts the alignments to BED format, and then counts how many hit each window in our 1kb genomic window file. 

Ok, now we've summarized coverage in units of reads per 1kb window. The file looks like this:

```bash
chr20	29399000	29400000	59	60	120
chr20	29400000	29401000	56.99429224	60	1752
chr20	29401000	29402000	58.54471931	60	2102
chr20	29402000	29403000	59.48537005	60	1743
```


Columns 1-3 are standard BED columns. Columns 4-6 are the mean and median of the mapping quality scores, and the count of reads mapping to the window for all three samples. A coarse way of looking at the distribution across windows is below:


```bash
module load htslib
tabix coverage_1kb.bed.gz chr20:29400000-34400000 | cut -f 6 | sort -g | awk '{x=100}{$1=int($1/x)*x}{print $1}' | uniq -c
```

This code uses `tabix` to extract only the windows we grabbed data for, cuts out the 6th column of the file (the count of reads), sorts it by numeric value, uses `awk` to round to the nearest 100, and then prints only uniqe values and their counts. 

You can see that the modal value is 1700 reads per window, and most values are clustered around that. We can think of the mode as being the expectated coverage for well-behaved parts of the genome. Some windows, however, have much lower coverage, and a handful have extreme deviations, up to 40x the modal value (73000 reads per window!). That amounts to a combined coverage of 10,000x. 

You could load this file into R to make a plot (not part of this tutorial!) that could help you see a little better what's going on:

![windowed coverage](/img/coverage_plot.png)

We can see that some regions have excesses or deficiencies of coverage, and that these windows tend to cluster. In some cases, it is worth trying to sort out what's going on with coverage abberations using specific tools, as they may result from sample-specific structural variants (structural variation can be biologically [__very important__](https://www.nature.com/articles/nature15394)), but for short variant calling, we want to ignore them. 

Unfortunately, there's no discrete boundary between normal variation in coverage of single-copy genome sequence and problematic regions, so here we'll just keep windows that are less than 0.5 times the modal value and greater than 1.5 times the modal value. 

We can select those windows and merge them together like this:

```bash
zcat ../coverage_stats/coverage_1kb.bed.gz | \
awk '$6 < 850 || $6 > 2250' | \
bedtools merge | \
bgzip >../coverage_stats/coverage_outliers.bed.gz 
```

Note the last step here. You can a pipe text stream to bgzip so that it is immediately compressed. 

Using these thresholds we will exclude 672kb of our 5mb region. 

It's worth mentioning that `freebayes` can accept a copy number map, in BED format, that gives the copy number per sample, per region for the whole genome. If you identified CNVs of interest and wanted to call genotypes within them (or on the sex chromosomes, or mitochondrion), you could use that information in principle, rather than excluding those regions. 

___
scripts:
- [scripts/Part4a_coverage.sh](scripts/Part4a_coverage.sh)

## Call variants

Now that we've decided which regions to exclude from analysis, we can move forward with variant calling. Here we want to construct a pipeline that allows us to filter out the reads from offending regions and pass the rest to freebayes. There are several ways to do this, but we'll just look at one. 

```bash
# set a variable for the reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

OUTLIERWINDOWS=../coverage_stats/coverage_outliers.bed.gz

# note that bamtools region specification uses ".." instead of "-"
bamtools merge -list bam.list -region chr20:29400000..34400000 | \
bamtools filter -in stdin -mapQuality ">30" -isProperPair true | \
bedtools intersect -v -a stdin -b $OUTLIERWINDOWS | \
freebayes -f $GEN --stdin | \
bgzip -c >../variants_freebayes/chinesetrio_fb.vcf.gz
```

In this pipe we merge all three sample bams into a single stream using `bamtools merge`. We send that stream to `bamtools filter` to exclude reads with poor mapping quality and those that aren't properly paired. The reads are then passed to `bedtools intersect` with the `-v` flag set, and our set of outlier windows. In this part, any reads that overlap the windows we flagged above will be excluded. This filtered stream of reads is then passed to `freebayes` to actually call variants, and then the resulting VCF file is passed to `bgzip` to be compressed. 

So in this approach, instead of telling `freebayes` which specific regions to evaluate, we are instead letting it run over the whole region but excising the unwanted data from the stream. 

If we compare this approach to simply running `freebayes` over the whole region and filtering afterward, for example:

```freebayes
freebayes -f $GEN -0 -L bam.list -r chr20:29400000-34400000 | bgzip -c >out.vcf.gz
```

We see that excluding the regions with high and low coverage yields a run time of about 29 minutes for this region, while simply running it on the entire region results in a run time of 45 minutes. Most of that 50% increase in run time is attributable to the very small number of windows with very high coverage. 

___
scripts:
- [scripts/Part4b_freebayes.sh](scripts/Part4b_freebayes.sh)

