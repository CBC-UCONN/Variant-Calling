# Variant discovery tutorials

This repository contains a variant detection tutorial for UConn CBC workshops. 

## Introduction

This repository an introduction to the basics of variant calling from high-throughput, short-read sequencing data. While some limited conceptual ground will be covered in the tutorial, if you are working through it independently, it will be much more helpful if you understand the motivation for the steps in advance. A useful (if dated) review of the underlying concepts is [Nielsen et al. 2011](https://www.nature.com/articles/nrg2986) in Nature Reviews Genetics. 

Most steps have associated bash scripts tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) scheduler. These can be modified to run interactively or with another job scheduler.  

Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should submit each script to the scheduler as `sbatch scriptname.sh` after modifying it as appropriate.  

Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this tutorial, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, you can get one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

__Structure:__

1. [ Stepwise QC, alignment, post-alignment processing ](/Part1_qc_alignment.md)

2. [ Variant Calling: bcftools ](/Part2_bcftools.md)

3. [ Part 1, but a piped example ](Part3_pipedalignment.md)

4. 
	a. [ Variant calling: Freebayes ](Part4a_freebayes.md)

	b. [ Variant calling: GATK, joint calling using gvcf ](Part4b_gatk.md)

5. [ Filtering and comparing variant sets ](Part5_filtering_comparing.md) 

6. [ Variant annotation ](Part6_annotation.md)

__Data:__

In this tutorial we will use human __whole genome shotgun sequence data__ from the [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project. The data are Illumina paired-end sequence data from a chinese trio consisting of a mother, father and son. For the son, the data are 2x250bp reads at 100x coverage. For the parents, the data are 2x150bp reads at 50x coverage. To make things run quickly, we'll only analyze a 5mb region of chromosome 20: chr20:29400000-34400000

This 5mb region of the genome is somewhat arbitrary, but it is near the centromere and has a few regions with mapping problems that can help illustrate technical issues that might arise during variant calling. 

For the reference genome, we'll use GRCh38. The specific version we'll use is [recommended by Heng Li, author of bwa](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) for variant calling. 

_Source:_
- [GiaB](https://www.nist.gov/programs-projects/genome-bottle)
- [GiaB FTP site](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/)
- Heng Li [recommends](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) to use this version of the human genome for variant calling:
	- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

__Required software tools:__

_quality control_  
- [ FastQC ](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ sickle ](https://github.com/najoshi/sickle)  

_alignment_  
- [ bwa ](http://bio-bwa.sourceforge.net/)

_sequence alignment manipulation and visualization_ 
- [ samtools ](http://www.htslib.org/doc/samtools.html)
- [ picard ](https://broadinstitute.github.io/picard/)
- [ samblaster ](https://github.com/GregoryFaust/samblaster)
- [ bamtools ](https://github.com/pezmaster31/bamtools)  
- [ igv ](https://software.broadinstitute.org/software/igv/)

_variant calling_  
- [ bcftools ](http://www.htslib.org/doc/bcftools.html)
- [ freebayes ](https://github.com/ekg/freebayes)
- [ GATK ](https://software.broadinstitute.org/gatk/)  

_other utilities_  
- [ bgzip ](http://www.htslib.org/doc/bgzip.html)
- [ tabix ](http://www.htslib.org/doc/tabix.html)
- [ bedtools ](https://bedtools.readthedocs.io/en/latest/)
- [ vcflib ](https://github.com/vcflib/vcflib)
- [ vt ](https://github.com/atks/vt)
