# Variant annotation

## Introduction

This section of the tutorial deals with annotating variants.  It assumes you have completed Parts 2, 4a/b and 5 so that you have 3 VCF files from each of the three variant calling approaches in the directory structure established in the tutorial. 

Steps here will use the following software packages:

- [ SnpEff ](http://snpeff.sourceforge.net/SnpEff.html)
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Update your working directory ](#Update-your-working-directory)  
3.    [ Predicting the functional effects of variants ](#Predicting-the-functional-effects-of-variants)

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

One of the most widely used tool for predicting functional effects of variants given a VCF and genome annotation is [`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), which is what we'll use here (but see also [VEP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4), and [ANNOVAR](https://academic.oup.com/nar/article/38/16/e164/1749458)). `SnpEff's` predictions stem entirely from the annotation of the primary DNA sequence, and so they do not account for the effects of secondary or tertiary structure of amino acid sequences, location with respect to important motifs within genes, or sequence conservation among species. For an approach that makes predictions based on amino acid sequence conservation across homologs, see [PROVEAN](http://provean.jcvi.org/index.php)

___

To run, `SnpEff` requires a database generated from a reference genome and annotation. There are many pre-built databases, which can be downloaded. The list can be viewed by typing:

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
You can refer to specific pre-built databases by the names in the first column, e.g. hg38. If the species you work on does not have a pre-built database, you can [build one](http://snpeff.sourceforge.net/SnpEff_manual.html#databases). 

Distributions of `SnpEff` do not contain all the pre-built databases, so if you're using one when you run it for the first time, it's likely it will need to download it. If you have your own copy of SnpEff, it will download it without a hitch. If you're using the `module` version on the Xanadu cluster, you will need to specify a full path to a data directory of your choosing with the option `-dataDir`, because you don't have write access to `SnpEff's` directory. 

That said, it's easy to run the program like this:

```bash
module load htslib/1.7
module load snpEff/4.3q

# freebayes filtered vcf file
VCF=../filtered_vcfs/fb_filter.vcf.gz

java -Xmx8G -jar $SNPEFF eff -dataDir $(pwd)/snpeff_data hg38 $VCF | bgzip -c > fb_filter.ann.vcf.gz
	
```

This will annotate our filtered `freebayes` output. `SnpEff` is very fast, so the longest part of this analysis will actually be downloading the 0.5G database. 

The output will be the original VCF file with the annotation added to the INFO field. Let's look at a couple that have impact categorized as "HIGH":

```bash
bcftools view -H -f PASS fb_filter.ann.vcf.gz | grep HIGH
```

With the resulting two VCF records:


```bash
chr20	31389027	.	A	C	1060.99	PASS	AB=0.483516;ABP=3.22506;AC=1;AF=0.166667;AN=6;AO=45;CIGAR=1X;DP=225;DPB=225;DPRA=2.11628;EPP=13.8677;EPPR=3.11948;GTI=0;LEN=1;MEANALT=1.5;MQM=60;MQMR=59.9553;NS=3;NUMALT=1;ODDS=41.7753;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=1577;QR=6054;RO=179;RPL=24;RPP=3.44459;RPPR=3.11948;RPR=21;RUN=1;SAF=24;SAP=3.44459;SAR=21;SRF=100;SRP=8.36013;SRR=79;TYPE=snp;ANN=C|stop_lost|HIGH|DEFB119|DEFB119|transcript|NM_153323.4|protein_coding|2/2|c.265T>G|p.Ter89Gluext*?|431/490|265/267|89/88||,C|stop_lost|HIGH|DEFB119|DEFB119|transcript|NM_001271209.1|protein_coding|2/2|c.262T>G|p.Ter88Gluext*?|428/487|262/264|88/87||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NM_153289.3|protein_coding|1/1|c.61+1396T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073151.1|pseudogene|1/3|n.228-668T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073152.1|pseudogene|1/2|n.228-668T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073153.1|pseudogene|1/2|n.227+1396T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_126440.1|pseudogene|1/1|n.227+1396T>G||||||	GT:DP:AD:RO:QR:AO:QA:GL	0/1:91:46,44:46:1452:44:1567:-114.109,0,-103.744	0/0:91:90,1:90:3078:1:10:0,-26.3937,-276.181	0/0:43:43,0:43:1524:0:0:0,-12.9443,-137.401
chr20	33217596	.	ACC	AC	1436.76	PASS	AB=0.504202;ABP=3.02855;AC=2;AF=0.333333;AN=6;AO=60;CIGAR=1M1D1M;DP=205;DPB=185;DPRA=0.69186;EPP=3.58936;EPPR=3.02549;GTI=0;LEN=1;MEANALT=1.5;MQM=60;MQMR=60;NS=3;NUMALT=1;ODDS=75.9105;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=2084;QR=4906;RO=143;RPL=26;RPP=5.32654;RPPR=4.8477;RPR=34;RUN=1;SAF=26;SAP=5.32654;SAR=34;SRF=80;SRP=7.3988;SRR=63;TYPE=del;ANN=AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_178466.4|protein_coding|1/7|c.62delC|p.Pro21fs|289/1181|62/765|21/254||,AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_001042439.2|protein_coding|1/6|c.62delC|p.Pro21fs|289/1073|62/657|21/218||;LOF=(BPIFA3|BPIFA3|2|1.00)	GT:DP:AD:RO:QR:AO:QA:GL	0/1:77:35,41:35:1207:41:1392:-102.499,0,-85.8728	0/0:86:85,0:85:2890:0:0:0,-25.5875,-259.891	0/1:42:23,19:23:809:19:692:-49.9554,0,-60.4688

```

Deep in the INFO field are the annotations, tagged with ANN:

```bash
ANN=C|stop_lost|HIGH|DEFB119|DEFB119|transcript|NM_153323.4|protein_coding|2/2|c.265T>G|p.Ter89Gluext*?|431/490|265/267|89/88||,C|stop_lost|HIGH|DEFB119|DEFB119|transcript|NM_001271209.1|protein_coding|2/2|c.262T>G|p.Ter88Gluext*?|428/487|262/264|88/87||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NM_153289.3|protein_coding|1/1|c.61+1396T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073151.1|pseudogene|1/3|n.228-668T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073152.1|pseudogene|1/2|n.228-668T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_073153.1|pseudogene|1/2|n.227+1396T>G||||||,C|intron_variant|MODIFIER|DEFB119|DEFB119|transcript|NR_126440.1|pseudogene|1/1|n.227+1396T>G||||||
ANN=AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_178466.4|protein_coding|1/7|c.62delC|p.Pro21fs|289/1181|62/765|21/254||,AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_001042439.2|protein_coding|1/6|c.62delC|p.Pro21fs|289/1073|62/657|21/218||
```

A VCF record can have more than one annotation, so each ANN field is a comma-separated list. The second record (position 33217596) has two annotations, corresponding to effects on two different transcripts of the gene BPIFA3:

```bash
AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_178466.4|protein_coding|1/7|c.62delC|p.Pro21fs|289/1181|62/765|21/254||
AC|frameshift_variant|HIGH|BPIFA3|BPIFA3|transcript|NM_001042439.2|protein_coding|1/6|c.62delC|p.Pro21fs|289/1073|62/657|21/218||
```

Each annotation is formatted as 16 subfields, separated by "|" and defined in the [manual](http://snpeff.sourceforge.net/SnpEff_manual.html). The first three are most salient: the allele, the effect (here, a frameshift for both annotations), and a predicted impact (here "HIGH"). The rest of the subfields describe the details of the feature and how the variant alters it. The effect and the feature type (subfield 6) are part of a controlled vocabulary, or ontology, called [the Sequence Ontology](http://www.sequenceontology.org/). 

`SnpEff` also produces the files `snpEff_summary.html` and `snpEff_genes.txt`. These summarize the set of variants that have been annotated. Download `snpEff_summary.html` to your local machine and view it in a web browser to see variant annotations broken down in various ways. 

Once you have annotated your VCF, you _can_ use `grep` to parse out variants with annotations you're interested in, as we did above, but [SnpSift](http://snpeff.sourceforge.net/SnpSift.html), which we won't cover here, may be a better approach. 

___
scripts:
- []()


## Connecting variants with existing databases



scripts:
- []()




bcftools annotate does simple matching. does not account for complex alleles or incompatible representations. 
`
--dbsnp in haplotypecaller can accept a vcf file and annotate snps

