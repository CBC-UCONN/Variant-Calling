# Variant annotation

## Introduction

This section of the tutorial deals with annotating variants.  It assumes you have completed Parts 2, 4a/b and 5 so that you have 3 VCF files from each of the three variant calling approaches in the directory structure established in the tutorial. 

Steps here will use the following software packages:

- [ SnpEff ](http://snpeff.sourceforge.net/SnpEff.html)
- []()

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  

## Motivation

Once you have identified a set of variants of interest, you may want to annotate them with information derived from outside resources. Two kinds of annotations are commonly applied to VCF files. The first annotation of functional prediction. You may want to know if any of the variants you have identified are likely to impact how the genome works. Do any occur in coding or promoter regions? Do any cause amino acid substitutions or premature stop codons? We do this by comparison with a functional annotation of the genome. Second, you may want to learn anything that is already known about variants you have found. Have alleles you observe previously been associated with phenotypic variation or disease? Here we will use both approaches. We do this by connecting them to databases of existing information (such as dbSNP). 

## Update your working directory

First we'll create a new directory to hold the results we'll generate. Make sure you're in the directory vc_workshop and type:

```bash
mkdir -p annotated_vcfs
cd annotated_vcfs
```

## Predicting the functional effects of variants

How does functional prediction work? Up until this point our input data has consisted of 1) a reference genome and 2) fragments of whole genome shotgun sequence from samples of interest. We have used that data to identify sequence variation among the sampled genomes. In order to understand whether that sequence variation has a biological impact, we could do experiments (not feasible for millions of variants in eukaryotic genomes!), or we could use sequence context to make some predictions. In this case, the sequence context is a **genome annotation**. 

A genome annotation can give us a wide array of information about genomic features of interest, but usually incorporates at least the location and coding frame of exons within genes, possibly for multiple transcripts. It may also identify non-coding RNAs, promoters, enhancers, conserved elements or repetitive elements. 

For a given genome annotation, variants can be categorized based on their predicted biological impact. As an aside, it's important to recognize that the quality of functional prediction is conditional on the quality and completeness of the genome annotation. 

One of the most widely used tool for predicting functional effects of variants given a VCF and genome annotation is [`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), which is what we'll use here (but see also [VEP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4), and [ANNOVAR](https://academic.oup.com/nar/article/38/16/e164/1749458)). `SnpEff's` predictions stem entirely from the annotation of the primary DNA sequence, and so they do not account for the effects of secondary or tertiary structure of amino acid sequences, location with respect to important motifs within genes, or sequence conservation among species. 

To run, `SnpEff` requires a database generated from a reference genome and annotation. There are many prebuilt databases, which can be downloaded. The list can be viewed by typing:

```bash
module load snpEff/4.3q
java -jar snpEff.jar databases
```
Which yields a list formatted like this:

```bash
hg19                                                        	Homo_sapiens (USCS)                                         	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip
hg19kg                                                      	Homo_sapiens (UCSC KnownGenes)                              	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19kg.zip
hg38                                                        	Homo_sapiens (USCS)                                         	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38.zip
hg38kg                                                      	Homo_sapiens (UCSC KnownGenes)                              	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38kg.zip
```
You can refer to specific pre-built databases by the names 


takes in a genome annotation, a reference genome, and a VCF file, and outputs a VCF file with annotations added to the INFO field, and tables of summaries of the predicted biological impact of the variants. 





___
scripts:
- []()


## Connecting variants with existing databases

The final thing we'll cover in this tutorial is variant annotation. In this context, annotation means two things:
1. Connect the variants we've discovered to knowledge about them in existing databases. 
2. Make some predictions about the functional impact of the variants based on our genome annotation. 

scripts:
- []()




bcftools annotate does simple matching. does not account for complex alleles or incompatible representations. 

--dbsnp in haplotypecaller can accept a vcf file and annotate snps

