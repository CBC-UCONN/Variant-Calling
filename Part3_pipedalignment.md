# Alignment, post-alignment processing using a pipe

## Introduction

This section of the tutorial revisits the alignment and post-processing steps from [Part 1](/Part1_qc_alignment.md), but introduces the concept of using unix pipes to streamline the workflow. 

We'll assume we have already downloaded and QC'ed fastq files as in steps 1-5 from Part 1 and that they exist in the directory structure first created there. 

Steps here will use the following software packages:

- [ samtools ](http://www.htslib.org/doc/samtools.html)
- [ bwa ](http://bio-bwa.sourceforge.net/)
- [ samblaster ](https://github.com/GregoryFaust/samblaster)

Each major step below has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 

## Contents
  
1.    [ Motivation ](#Motivation)
1.    [ Update your working directory ](#Update-your-working-directory)  
2.    [ Piped alignment and post-processing ](#Piped-alignment-and-post-processing)
3.    [ Index alignment files ](#Index-alignment-files)

## Motivation

In [Part 1](/Part1_qc_alignment.md) of the tutorial, we downloaded data, QC'ed it, aligned it to a reference genome, and did post-alignment processing. Each step was conducted using a discrete script. Each script required that all the data be read, modified, and then re-written. This creates two problems. First, for large datasets, all of this reading and writing takes time. Second, it generates many copies of the data. If our first approach had been applied to a 300 gigbyte dataset, by the time we finished we'd be using somewhere in the neighborhood of 1.8 terabytes of disk space. 

In this section of the tutorial we'll see how we can streamline this approach. To do this we'll make use of a core feature of unix-like operating systems: the pipe. 

The pipe, `|`, allows you to redirect the output of one command to the input of another. In unix-speak, you would say that you are [redirecting the "standard output" stream of command 1 to the "standard input" of command 2](https://en.wikipedia.org/wiki/Standard_streams). 

In this way, you can chain together many specialized programs to achieve some goal without having to repeatedly read and write files. This mitigates both of the problems described above. 

Many programs in bioinformatics are designed with pipes in mind, including core utilities such as `samtools` and the variant caller `Freebayes` that we'll use in Part 4. 


## Update your working directory

First we'll make a new directory to hold our second set of alignments. 

```bash
 mkdir align_pipe
 ```

## Piped alignment and post-processing

Our goal here in using the pipe is to read and write the data as few times as possible, within reason. So each time we read the data, we want to accomplish as many tasks with it as we can before we have to write it again. Some tasks are better suited to this than others. For example, alignment occurs on each read pair independently, so as we read through our data, the aligner can do its work and pass the alignment off to the next step. Compression works similarly well, as the compression algorithms we use don't need to see the entire file at once to do their work. Sorting and marking duplicates is a little more complicated, and the software we use for those will write temporary files if there is too much data to hold in memory, but piping still makes things more efficient. 

Here is some pseudocode to give an idea of what what we're aiming at here:


```bash
align | mark duplicates | compress | sort >sample.bam.
```
In our first step we'll align the data. In the second we pass it to the duplicate marking software. In the third we compress it. In the fourth we sort it, and finally write it to a file. 

For example:

```bash
# set a variable 'GEN' that gives the location and base name of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38

# execute the pipe:
bwa mem -t 4 -R '@RG\tID:son\tSM:son' $GEN ../rawdata/son.1.fq ../rawdata/son.2.fq | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/son - >son.bam
```

Here we run `bwa mem` the same way as in Part 1. We pass the alignments to `samblaster`, which marks duplicates (we used `picard tools` previously). `samblaster` doesn't need any options specified in this case. We then pass the aligned, duplicate-marked sequences to samtools to first compress, then sort them. For samtools the `-` in each command indicates that the data should be read from the standard input. For `samtools sort`, the `-T` indicates where temporary files should be written. Depending on the architecture of a given computer cluster, it may be better to write temporary files to a special disk, often titled /scratch, which may have faster read and write speeds than your home directory. Finally, we write the alignments to a bam file. 

Instead of reading and writing our sequences 5 times, resulting in 6 copies of our sequence data, we now only have two copies (or 3 if we have quality trimmed), the original fastq files, and processed bam files. 

If we had multiple bam files for each sample, we could then merge them into a single file. this is convenient, although not strictly necessary: both GATK and Freebayes identify which reads belong to which samples by the read group tags, and not the bam file from which they originated. 

___
scripts:
- [scripts/Part3a_alignment.sh](scripts/Part3a_alignment.sh)

## Index alignment files

As in Part 1, the last step in preparing the reads is to index the bam files. This needs to be done to enable fast access to reads overlapping a given position of the genome. Without the index, if you wanted to access reads at the beginning of chromosome 20, you'd need to read through chromosomes 1-19 until you got there. With many samples or deep coverage, this would be a big problem. The bam index is a map to the bam file that lets you skip around quickly. We'll use samtools to do this. For example:   

```bash
samtools index ../align_pipe/*.bam
```

scripts:
- [scripts/Part3b_indexbams.sh](scripts/Part3b_indexbams.sh)
