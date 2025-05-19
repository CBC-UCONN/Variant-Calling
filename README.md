# Variant discovery tutorials

This repository contains a variant detection tutorial for UConn CBC workshops. 

## Introduction

This repository an introduction to the basics of variant calling from high-throughput, short-read sequencing data. While some limited conceptual ground will be covered in the tutorial, if you are working through it independently, it will be much more helpful if you understand the motivation for the steps in advance. A useful (if dated) review of the underlying concepts is [Nielsen et al. 2011](https://www.nature.com/articles/nrg2986) in Nature Reviews Genetics. 

The tutorial assumes that you are working on a high performance computing cluster in a Linux environment. Most steps have associated scripts written in [bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) and specifically tailored to the UConn CBC Xanadu cluster. They have headers for the [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler and use [environment modules](http://modules.sourceforge.net/) to load software. 

These scripts should be easily modifiable to run without using SLURM or pre-installed software modules if you would like to try this out somewhere other than Xanadu. 

If you are working on Xanadu, commands should never be executed on the submit nodes of any HPC machine.  You should submit each script to the scheduler as `sbatch scriptname.sh` after modifying it as appropriate, or request an interactive session (e.g. `srun -c 2 --mem=25G --qos=general -p general --pty bash`)

Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this tutorial, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, you can get one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

__Data:__

In this tutorial we will use human __whole genome shotgun sequence data__ from the [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project. The data are Illumina paired-end sequence data from a chinese trio consisting of a mother, father and son. For the son, the data are 2x250bp reads at 100x coverage. For the parents, the data are 2x150bp reads at 50x coverage. To make things run quickly, we'll only analyze a 5mb region of chromosome 20: chr20:29400000-34400000

This 5mb region of the genome is somewhat arbitrary, but it is near the centromere and has a few regions with mapping problems that can help illustrate technical issues that might arise during variant calling. 


_Source:_
- [GiaB](https://www.nist.gov/programs-projects/genome-bottle)
- [GiaB FTP site](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/)

