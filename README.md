# Variant discovery tutorials

This repository contains re-worked variant detection tutorials for UConn CBC workshops. 

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

	c. Beyond variant calling: genotype likelihoods. 

5. [ Filtering and comparing variant sets ](Part5_filtering_comparing.md) 

6. [ Variant annotation ](Part6_annotation.md)

__Proposed data:__

[NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) chinese trio. 
Region chr20:29400000-34400000

This is an arbitrary 5mb region of the genome. It is:
- near the centromere. 
- has a couple regions with mapping problems. 
- takes only a few minutes to align. 
- 50-100x whole genome shotgun sequencing for each of 3 individuals. 
	- 2x250bp paired-end reads for the son
	- 2x150bp paired-end for the parents

For a reference, we'll use GRCh38. The specific version we'll use is [recommended by Heng Li, author of bwa](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) for variant calling. 

_Source:_
- [GiaB](https://www.nist.gov/programs-projects/genome-bottle)
- ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
- [Heng Li's recommendation](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
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
