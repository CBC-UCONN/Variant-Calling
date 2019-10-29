# Variant annotation

## Introduction

This section of the tutorial deals with annotating variants.  It assumes you have completed Parts 2, 4a/b and 5 so that you have 3 VCF files from each of the three variant calling approaches in the directory structure established in the tutorial. 

Steps here will use the following software packages:

- []()
- []()

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  

## Motivation

Once you have identified a set of variants of interest, you may want to annotate them with information derived from outside resources. Two kinds of annotations are commonly used. The first is functional prediction. Conditioned on a functional annotation of your reference genome, you may want to know if any of the variants you have identified are likely to impact the functioning of the genome. Do any occur in coding or promoter regions? Do any cause amino acid substitutions? Second, you may want to learn anything that is already known about variants you have found by connecting them to databases of existing information (such as dbSNP). Here we will use both approaches. 

## Update your working directory

First we'll create a new directory to hold the results we'll generate. Make sure you're in the directory vc_workshop and type:

```bash
mkdir -p annotated_vcfs
cd annotated_vcfs
```

## Predicting the functional effects of variants

Perhaps the most common 

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

