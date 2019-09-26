# Variant discovery using bcftools. 

## Introduction

This section of the tutorial introduces variant calling using the methods implemented in `bcftools`. 

We will assume that you have completed Part 1 and have QC'ed and processed sequence alignment files in bam format in the directory structure first created there. 

Steps here will use the following software packages:

- [ bcftools ](http://www.htslib.org/doc/bcftools.html)
- [ tabix ](http://www.htslib.org/doc/tabix.html)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Generate a pileup file ](#Generate-a-pileup-file)  
3.    [ Call variants ](#Call-variants)  
4.    [ The VCF format ](#The-VCF-format)  

## Motivation

In [Part 1](/Part1_qc_alignment.md), we conducted the first phase of the reference-mapping approach to variant discovery: we aligned short read sequencing data to a reference genome and prepared the alignment files for variant calling. In this part, we'll see how we evaluate the resulting sequence alignments for evidence of variation using `bcftools`, which applies a statistical model to evaluate the evidence for variation. 

Why do we apply a statistical model to discover variants? Why can't we just say "this observed sequence differs from the reference genome at this position, therefore we have found a variant"?

In short, because the sequence data we observe have been passed to us through several messy and error-prone laboratory and computational processes. We can break down this pipeline into several steps:

1. _Sampling the genome:_ A diploid individual has two copies of each chromosome. If this individual is heterozygous for a given site, the probability that we only observe the reference allele is 0.5^n, where n is the number of (non-duplicate) sequences. We still have a 3.125% chance of failing to observe the alternate allele given 5x coverage of the site. 
2. _Library preparation:_ During library preparation, polymerases can misincorporate bases. These bases will be read by the sequencer and likely assigned high base quality scores. If PCR is used during library preparation, random changes in the allele balance at heterozygous sites can be introduced, or existing random biases from the pool of extracted DNA can be amplified, exacerbating the problem in the previous step. 
3. _Sequencing:_ When we load our library of DNA fragments onto the sequencer to be analyzed, the machine may read bases incorrectly. 
4. _Reference mapping:_ After sequencing our library, we must map the sequences back to the reference genome. Most reference genomes are a) incomplete and b) have many copies of similar or identical sequences. These issues can mean that a DNA fragment mapped back to the reference genome may not actually have been derived from that region. Differences between it and the reference therefore don't represent the kind of polymorphism we're looking for. 
5. _Sequence alignment:_ Even if we have assigned a DNA fragment to the correct location in the genome, we still need to align it properly. If there are indels, this can become very difficult and there may be several equally likely alignments. 

So, all of this is to say that when we observe a set of sequences aligned to the reference genome, and there appear to be some bases that differ, we need a statistical model that can account for these complex sources of error to help us decide whether that variation is real and assign genotypes to individuals. 

## Update your working directory

First we'll make a new directory to hold variant calls. 

```bash
mkdir variants_bcftools
```

## Generate a pileup file

`bcftools` uses a two step procedure to call variants. First it generates a pileup file using `bcftools mpileup`. The pileup file summarizes the per base coverage at each site in the genome. Each row represents a genomic position, and each position that has sequencing coverage is present in the file. The second step, using `bcftools call` actually applies the statistical model to evaluate the evidence for variation represented in that summary and generate output in the variant call format (VCF). 

Because pileup files for whole genome sequencing contain a summary for every site in the genome for every sample, they can be very large. We also don't typically need to look at them directly, or do anything with them outside of th variant calling. So if you use this approach on your own data, you will usually simply pipe the output of `bcftools mpileup` directly to `bcftools call` (for further discussion of the use of pipes, see section 3). 

For demonstration purposes here, we'll keep the steps separate. 

Here is a call to `bcftools mpileup`:

```bash
# set reference genome location
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta

bcftools mpileup \
	-f $GEN \
	-b list.bam \
	-q 20 -Q 30 \
	-r chr20:29400000-34400000 >../variants_bcftools/chinesetrio.pileup
```

We give `bcftools` the reference genome with `-f`, a list of bam files with `-b`, tell it to exclude bases with quality lower than 30 and reads with mapping quality lower than 20 with `-q` and `-Q` and ask it to generate the pileup only for the region we're focusing on here with `-r`. 

___
scripts:
- [scripts/Part2a_mpileup.sh](scripts/Part2a_mpileup.sh)

## Call variants

The second step is to evaluate the evidence (summarized in the pileup file) that the sequence variation we observe is true biological variation, and not errors introduced during library preparation, sequencing, mapping and alignment. Here we use `bcftools call`. 

```bash
bcftools call -m -v -Oz -o chinesetrio.vcf.gz chinesetrio.pileup
```

We use the `-o` flag to indicate the output file name. The flag `-m` specifies one of two possible variant calling routines, `-v` says that only variable sites should be output, and `-Oz` indicates the output should be compressed in a version of the gzip compression format. 

In this case, we've told bcftools to output a compressed file. Other variant callers may not have that option. In those cases, it's a good idea to use compression to save space. `bgzip` is a commonly used utility for compressing tabular data in genomics. It has the advantage of producing files that can be indexed, generally using the companion program `tabix`. bcftools used the same algorithm as bgzip to compress our file, so we'll use `tabix` to index it so that we can access quickly access variants from any part of the genome. 

This indexing is unnecessary for our 5mb chunk of the human genome, but it is critical for VCF files containing millions of variants. 

___
scripts:
- [scripts/Part2b_variantcall.sh](scripts/Part2b_variantcall.sh)
- [scripts/Part2c_tabix.sh](scripts/Part2c_tabix.sh)

## The VCF format

Now we have a set of VCF formatted variants. Before going further, we should learn a little about what's in this file. There are no scripts for this section, so execute the code yourself on the command line. 

We generated a block-gzip compressed VCF file in the last step. We will inspect this file using `bcftools view` So make sure `bcftools` is loaded by entering `module load bcftools` (Alternatively `zcat` will print it to the screen, and we can read it using `less`). 

The first 3000 or so lines of the file are the header, which contains a lot of metadata on the file. Each header line starts with `##`. To view the header say:

```bash
bcftools view -h chinesetrio.vcf.gz
```
The `-h` flag prints only the header. 

The first few header lines give you basic information about the format and origin of the file:

```bash
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.6+htslib-1.6
##bcftoolsCommand=mpileup -f /UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta -b list.bam -q 20 -Q 30 -r chr20:29400000-34400000
##reference=file:///UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
```
The biggest chunk of header lines lists the sequences in the reference genome and their length. In this case there are 3366 sequences, but the vast majority of the genome is in the first 24: 22 autosomes + X + Y. 

Finally the last 25 or so lines of the header give the definitions for abbreviations used in the variant records to follow, for example:

```bash
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
```

The "DP" tag in the INFO field (see below) gives the raw read depth (across samples) for the variant. 

The final header line gives a set of field names:

```bash
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	dad	mom	son
```

After the header come individual variant records. Each record, taking up one line, represents a possible variant.

To view the first record, try:

```bash
bcftools view -H chinesetrio.vcf.gz | head -n 1
```

It should result in a line like this:

```
chr20	10000775	.	A	G	222	.	DP=269;VDB=0.0730423;SGB=51.8218;RPB=0.806901;MQB=1;MQSB=1;BQB=0.929014;MQ0F=0;ICB=0.3;HOB=0.125;AC=1;AN=4;DP4=105,74,27,22;MQ=60	GT:PL	0/0:0,255,255	0/1:255,0,255	./.:0,0,0
```

The first two columns are the reference sequence and position of the variant. Then comes the "ID" field, now empty ("." in every entry) which you can population by comparison to a database or other VCF file (see Part 5). Next are the reference and alternate alleles. If more than one alternate allele is observed, there will be a comma-separated list of alleles. Field 6 contains the [phred-scaled](https://en.wikipedia.org/wiki/Phred_quality_score) variant quality score (i.e. 10^( QUAL / -10) = probability of false call). Then comes the filter field, now empty, which can be populated by PASS, or a list of filtering criteria the variant has failed. Next comes the INFO field, which generally contains many summaries of the read data supporting the variant, identified by abbreviations defined in the header. After this is the FORMAT field, which gives the formatting of the individual genotypes in the columns that follow. In this case `GT:PL`, means that the genotype is given first, followed by a colon, and then phred-scaled genotype likelihoods. The next three columns are the genotypes for son, mom, and dad, respectively. As an example, for the variant above, the mother's genotype is given by:  

`0/1:255,0,255`  

"0/1" means her genotype is called as heterozygous (0 is the reference allele, 1 is the alternate). The numbers 255,0,255 are the phred-scaled genotype likelihoods for genotypes 0/0, 0/1, 1/1 respectively. These are scaled by the likelihood of the best genotype. Because that makes the likelihood of that genotype 1, and phred is log-scaled, the best likelihood is always 0. This also means you can think of the other genotype scores as ratios of that genotype's likelihood to the best likelihood. Thus for a likelihood of 255, the probability of the observed data being produced by the heterozygous genotype is 3x10^25 times the probability of it having been produced by either homozygous genotype. 

For lots more information on VCF, [here's a link](https://samtools.github.io/hts-specs/VCFv4.2.pdf) to the format specification. 


