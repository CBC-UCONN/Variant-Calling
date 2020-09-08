# Stepwise QC, alignment, post-alignment processing. #

## Introduction

This section of the tutorial deals with quality control, alignment to a reference genome, and post-processing of high-throughput short-read sequencing data in preparation for variant calling. In this section, steps are executed individually. For an example using pipes, see Part 3. 

Steps here will use the following software packages:

- [ samtools ](http://www.htslib.org/doc/samtools.html)
- [ picard ](https://broadinstitute.github.io/picard/)
- [ bwa ](http://bio-bwa.sourceforge.net/)
- [ FastQC ](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ sickle ](https://github.com/najoshi/sickle)
- [ bedtools ](https://bedtools.readthedocs.io/en/latest/)
- [ igv ](https://software.broadinstitute.org/software/igv/)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or in other contexts. 


## Contents
  
1.    [ Motivation ](#Motivation)
1.    [ Set up a working directory ](#Set-up-a-working-directory)  
2.    [ Prepare a reference genome ](#Prepare-a-reference-genome)
3.    [ Download data ](#Download-data)
4.    [ Assess read quality ](#Assess-read-quality)
5.    [ Quality trim ](#Quality-trim)
6.    [ Align and compress ](#Align-and-compress)
7.    [ Sort reads by genome position ](#Sort-reads-by-genome-position)
8.    [ Mark duplicates ](#Mark-duplicates)
9.    [ Index alignment files ](#Index-alignment-files)
10.   [ Exploring SAM files ](#Exploring-SAM-files)

## Motivation

Reference mapping is the dominant paradigm for discovering variants using high throughput short read sequencing data. The general idea is to obtain short read sequences from either the whole genome or a set of target regions (e.g. through exome capture or RAD-seq), map those sequences back to a (hopefully) highly contiguous reference genome, and use a statistical model to discover variants and genotype individual samples. 

Reference mapping contrasts with other approaches, such as de novo genome assembly with whole genome alignment, which is currently far more expensive and computationally demanding. 

This section of the tutorial covers the first phase of the reference mapping approach: data QC, mapping, and post-processing of alignment files. 

## Set up a working directory ##

We will begin the tutorial by setting up a working directory to organize the files we'll generate. From your home directory, enter the code below on the command line to create a set of directories. 

```bash
mkdir -p vc_workshop/rawdata vc_workshop/fastqc vc_workshop/align_stepwise vc_workshop/scripts
cd vc_workshop
```


## Prepare a reference genome

The first step is the prepare the reference genome. Most software packages that align short-read sequencing data to, or otherwise manipulate a reference genome require that genome to be indexed in some way. We will generate indexes using both `bwa` and `samtools`. For the workshop, a pre-indexed human genome is provided, but the code below shows how you can download and index a human genome yourself:

```bash
# make bgzip available
module load htslib
# download a human genome, version as recommended by Heng Li: https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# compress the genome using bgzip
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | bgzip >GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# set a variable 'GEN' that gives the location of the reference genome:
GEN=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
bwa index $GEN
samtools faidx $GEN
```

If you are working on this tutorial somewhere other than UConn's xanadu cluster, you can make a directory called "reference", place this genome there, and edit following scripts to point at that version of the genome. 

## Download data ##

For all following steps of the workshop, we'll use data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project, hosted by the [National Institute of Standards and Technology](https://www.nist.gov). The data were collected as part of an effort to create reference standards to compare genomic workflows. Large quantities of sequencing data have been generated on multiple platforms (Illumina, PacBio, etc) for seven individuals. 

We're going to use data from a trio (mother, father, and son) of Chinese ancestry. The data consist of 250bp paired end reads sequenced on an Illumina HiSeq 2500. To ensure the analyses run quickly, we'll only use data from 5 megabases of chromosome 20 (position 29,400,000 to 34,400,000). The expected sequencing coverage is 100x for the mother and father, and 45x for the son. 

More information about the data can be found at the links below:  
- https://www.nist.gov/programs-projects/genome-bottle
- ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
___

To download the data, we'll use `samtools`. The data have already been aligned to a reference genome, and the resulting BAM file has been compressed and indexed. This will allow us to get reads only from the region we're interested in. Conveniently, `samtools` can read BAM files from an ftp server, provided the index is present, so we won't need to download the whole dataset. We'll then convert the data back to the unaligned fastq format using `bedtools` so we can continue with the tutorial.

We'll accomplish this with a unix pipeline, where the symbol `|` indicates that the output of the command to the left should be redirected as input to the command to the right. This will be discussed in more detail in Part 3 of the tutorial. 

The commands to be chained together are as follows:

- `samtools view -bh <file> <region>` outputs the given `region` from the `file`, includes the header, and outputs as a compressed BAM file.  
- `samtools sort -n -` sorts the reads by name so that read pairs will be found together in the file. The `-` indicates the data should be read from the pipe.    
- `bedtools bamtofastq -i /dev/stdin/ -fq <forward reads> -fq2 <reverse reads>` converts bam format back to fastq. In this case `/dev/stdin/` indicates the data should be read from the pipe. 

The result will be two fastq files, one containing the forward reads and one containing reverse reads with members of each pair on the same lines of their corresponding files. 

Below is an example of this code put together to download data for the son:

```bash
# set a variable 'SON' that gives the location of the file containing the data from the son:
SON='ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_45x_bams_vcfs_PerFlowCell/150424_HG005_Homogeneity_02_FCA-22108087/150424_HG005_Homogeneity_FCA_Combined-23168145/150424-HG005-Homogeneity-FCA-Combined_S1.bam'
# download the data, sort it, reformat to fastq
samtools view -uh $SON chr20:29400000-34400000 | samtools sort -n - | bedtools bamtofastq -i /dev/stdin/ -fq son.1.fq -fq2 son.2.fq
```
___
scripts:
- [scripts/Part1a_datadownload.sh](scripts/Part1a_datadownload.sh)

### Inspecting fastq files ###

This set of commands will allow you to open the FASTQ file and inspect the data.  Each read is represented by 4 lines.  Details on the FASTQ file format can be found <a href="https://bio Informatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>.   
 

Change the directory to the `rawdata/` folder using:   
```bash
cd rawdata/  
```  

To inspect the first few lines in the FASTQ file we can use the `head` command, as follows which will print the first few lines into the terminal window.     

```bash
head son.1.fq
```
A single fastq sequence entry looks like this:
```
@H2YVCBCXX:1:1:383:0/1
GCCCTGCCACACTTAGTGTATTGGCTTTTGTCTTCAGGGTTTTCACCTCGTGTTTTTGATATGGCTACCACGACTCCTAGACTACCTCATTACACAGCATTTACAGGCCAGGATAGGCTCCTCTTCATGTCATTATATGACAATAATTTTCATGTCCTTTTTAAAGCATAAAATAAAAATTCCCACACACCTCCAATAAAACTTTCTGTTCTGTCCCATAAGCTAGAACTGCATCTCTTGCCCATTCCT
+
DDDDDIIIIIHIHIIIEHHIIIIIIIGHHIIIIIIIIIGGHHIIHHIIIIIIHIGIHHHIHIFIHIHHHHIGHGHIIHHIHHIHIIIIHIIIIIIGIIEFHDCHHHIIIHIHHIIIHIHIIIIIIIHIGIIGHHIIIGEHHGGIGGHIECEHHIIIIIIHIGHH@H@GHHHEHHHHHIIF@@F?CHIH/CGGEHH@G/FHHHHHHFHHHHE@.DCHHHE..:GHE@.8.ABEHHEAHEG.BG.9B@BG.
```

The first line is the sequence name. The second line is the DNA sequence. The third line, which always begins with a "+" can contain other optional information, but usually does not. The fourth line encodes the base quality scores in ASCII characters. 

The command below will let you inspect more of the file.
```bash
less son.1.fq
```   

The command below will count the number of reads in the file. 
```bash
grep -c "^+" son.1.fq
```   

It's worth noting that while we are using uncompressed fastq files in this example, all of the programs we are working with will accept files compressed using gzip (generally denoted with the file extension '.gz'), and it's good practice to keep fastq files compressed so they use less storage space. 


## Assess read quality ##

In order to evaluate the general quality of reads in the file we will use the `FastQC` package.  The command can take multiple files as input and outputs a quality report in html format. 

Our script runs fastqc as follows:

```bash
fastqc -t 4 -o ../fastqc ../rawdata/*fq
```

The `*fq` indicates that fastqc should run on all files ending in "fq" in directory `rawdata/`. 

Once the files are generated you'll have to transfer them to your local computer to open them and examine the results. To copy the file from the Xanadu cluster please use the `transfer.cam.uchc.edu` node.  

You can use an ftp program, or the unix utility `scp` as:
```bash
scp user_name@transfer.cam.uchc.edu:/FULL_PATH_to_FILES/*.html . 
```

Again, `*html` will indicate that `scp` should copy all files ending in ".html". Once the files are downloaded, you can open them in a web browser. 

INSERT FIGURES HERE? OR LEAVE THEM TO WORKSHOP? 
___
scripts: 
- [scripts/Part1b1_fastqc.sh](scripts/Part1b_fastqc.sh)

## Quality trim ##

Current variant callers account for uncertainties in mapping (conditional on the quality of the reference genome) and in base calling, so quality trimming is not generally necessary for this application (the worrisome sources of error in variant calling are ["unknown unknowns"](https://en.wikipedia.org/wiki/There_are_known_knowns), like the incompleteness of the reference genome, or systematic error arising from library prep). However, if you have a dataset plagued by adapter contamination or poor quality reads, you may want to try trimming to salvage it and/or remove some of the noise. 

`sickle` is a commonly used tool for this task. 

We could trim the reads for the son as follows:

```bash
SEQ=son
sickle pe -t sanger \
	-l 100 \
	-f ../rawdata/$SEQ.1.fq \
	-r ../rawdata/$SEQ.2.fq \
	-o ../rawdata/$SEQ.trim.1.fq \
	-p ../rawdata/$SEQ.trim.2.fq \
	-s ../rawdata/$SEQ.trim.0.fq
```

This would discard any read trimmed shorter than 100bp, and if its pair was longer than 100bp, it would be placed in the file given by `-s`, which would be read as `../rawdata/son.trim.0.fq`. 

We can run our script, which also runs `FastQC` on the trimmed data, but as we will see, it will have little impact for this particular dataset. The next steps will use the untrimmed data. 
___
scripts:	
- [scripts/Part1b2_sickle_fastqc.sh](scripts/Part1b2_sickle_fastqc.sh)    


## Align and compress ##

Now that we have QC-ed our sequence data, it's time to align it to a reference genome. For that we'll use `bwa`, one of the most widely used short-read aligners. `bwa` implements several alignment methods, but `mem` is best for our application. We previously indexed our reference genome, so we're ready to go here. 

```bash
# set a variable 'GEN' that gives the location and base name of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38
# execute bwa mem
bwa mem -t 4 -R '@RG\tID:son\tSM:son' $GEN ../rawdata/son.1.fq ../rawdata/son.2.fq -o ../align_stepwise/son.sam
```
Here `-t 4` indicates that the program should use four processors, and we feed `bwa mem` both the location of the reference genome and both files of paired end reads.   

The `-R` flag specifies the read group information. Adding read group information is critical, though it does not have to be done at this stage (picard has a tool for adding read groups to alignment files). The read group can specify the source of the reads, including the library, sequencing machine, run, lane (for more see the SAM specification) but critically, for multisample variant calling it should include the sample ID. In this case we've specified the read group ID ("ID:") and the sample ID ("SM:") as the same thing. 

When simultaneously calling variants on many samples, variant callers do not track reads by their alignment file of origin, but using read group information. If this information is absent, all reads will be treated as if they came from a single sample, and a single genotype call at each variable site will result. Maintaining detailed read group information can also be helpful if some sequencing runs turn out to be problematic. In that case, even if all reads are pooled in a single bam file, there are tools you can use to filter out undesirable read groups on the fly.  

Finally, we'll compress the resulting alignment file:

```bash
samtools view -bhS ../rawdata/son.sam >../align_stepwise/son.bam
```

Because we're doing all steps individually, you will note that at this point (if we've done the quality trimming) we now have 4 copies of our sequence data. 
___
scripts:	
- [scripts/Part1c_align.sh](scripts/Part1c_align.sh)    
- [scripts/Part1d_compress.sh](scripts/Part1d_compress.sh)

## Sort reads by genome position ##

To call variants at a given position in the reference genome, we need to look at all the reads that overlap that position. In order to do this efficiently, we need to sort the reads in the alignment files by their positions in the reference genome. We'll use `picard tools` for this. For example:

```bash
IN=../align_stepwise/son.bam
OUT=../align_stepwise/son.sort.bam
java -jar $PICARD SortSam \
        INPUT=$IN \
        OUTPUT=$OUT \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True
```

scripts:	
- [scripts/Part1e_sort.sh](scripts/Part1e_sort.sh)

## Mark duplicates ##

Duplicate sequences are those which originate from the same molecule after extracting and shearing genomic DNA. There are two types: _optical_ and _polymerase chain reaction (PCR)_ duplicates. Optical duplicates are an error introduced by the sequencer. PCR duplicates are introduced by library prepartion protocols that use PCR. Duplicates cause 2 types of artifacts that mislead variant callers. 
- __First__, errors introduced by the polymerase can be propagated to multiple copies of a given fragment. Because these errors are actually part of a DNA sequence, they are likely to have high base qualities. If many sequences from the fragment containing the error are present, the variant caller can be deceived into identifying it as true biological variation. 
- __Second__, when variant callers call genotypes, they assume that heterozygous sites will have equal representation of both alleles in the sequence pool (as they should for germ-line mutations). Dramatically unbalanced coverage of an allele can be a signal that variation is spurious. Because of its exponential reproduction of fragments, PCR can randomly alter allele balance, or amplify small deviations in the initial sample, causing a variant caller to incorrectly call genotypes as homozygotes, or a truly variable site as invariant. 

For these reasons we need to exclude duplicate sequences from variant calling. They can be identified most easily from paired-end data as those sequences for which both reads have identical start sites. This may eliminate some sequences which are in fact derived from unique fragments in the original library, but if fragmentation is actually random, identical fragments should be rare. Once identified, duplicate sequences can be marked and ignored during variant calling (or other types of analyses) downstream. 

Here is some example code:

```bash
# assign input and output file names to variables
IN=../align_stepwise/son.sort.bam
OUT=../align_stepwise/son.mkdup.bam
# run MarkDuplicates
java -jar $PICARD MarkDuplicates \
        INPUT=$IN \
        OUTPUT=$OUT \
        METRICS_FILE=$IN.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=True
```

In this example, duplicate sequences remain in the file, but they are flagged as such. 

___

scripts:
- [scripts/Part1f_markduplicates.sh](scripts/Part1f_markduplicates.sh)

## Index alignment files ##

The last step in preparing the reads is to index the bam files. This needs to be done to enable fast access to reads overlapping a given position of the genome. Without the index, if you wanted to access reads at the beginning of chromosome 20, you'd need to read through chromosomes 1-19 until you got there. With many samples or deep coverage, this would be a big problem. The bam index is a map to the bam file that lets you skip around quickly. In this case, `picard` already generated the index during duplicate marking, but if you used a different approach, you'd need to use `samtools` to do this. For example:   

```bash
samtools index ../align_stepwise/*mkdup.bam
```

Now we have completed the initial QC, alignment and processing steps. At this point, you may have noticed that we have accumulated six copies of our data. Two copies of the fastq files, and four copies of the alignment files. This is a large and space-wasting mess. If we were working with many samples of high coverage human genomes, we would want to go and delete the intermediate alignment files and the trimmed fastqs, keeping only the original fastqs and the analysis-ready bams. Another approach, detailed in [Part 3](Part3_pipedalignment.md), would pipe many of these steps together and avoid creating some of the intermediate files to begin with. 

___

scripts:
- [scripts/Part1g_indexbams.sh](scripts/Part1g_indexbams.sh)


## Exploring SAM files ##

Now that we have processed alignment files we can learn about the sequence alignment format (SAM) and do some exploring. The ability to inspect and summarize various aspects of an alignment file can be especially helpful when something has gone wrong with sequencing or bioinformatic analysis. 

___


The best place to go for information about the SAM format is the [formal specification](https://samtools.github.io/hts-specs/SAMv1.pdf). But here we can discuss some of the main fields. Below are two entries from a SAM file:

```
HISEQ1:62:HB657ADXX:2:2106:14568:84364  1171    chr20   28557370        40      82M1D52M1I13M   =       28557055        -463    AGATGGTCTTGGTCCTTTTCCTGTTATGTGGAATCGTCAATTCAAATTTTAAAAGTGACTTTGAGATGTTTTTCATCTATTATCTAAAAATGTTGAAGGCGTTTTTAATTCTGCCTTCAACAGAGATCAACATACCTCTGATTATGAT    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDEDDDDDDDDEDDEDDDECDDDDDDDDDDCDEEEEEEEFFFFHHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJIJJJJJJJJJJJJJJJJJJIHGAJHHHHHFFFFFCCC    NM:i:8  MD:Z:5T2C73^T1T43G0G3C14        MC:Z:69S79M     AS:i:106        XS:i:111        RG:Z:dad        XA:Z:chr20,+30424698,63M1D85M,7;chr4,-189273762,82M1D66M,9;
HISEQ1:60:HB6D7ADXX:2:1116:1268:60355   83      chr20   28557379        40      4S126M1I17M     =       28557054        -468    GTCTTGGTCCTTTTCCTGTTATGTGGAATCGTCAATTCAAATTTTAAAAGTGACGTTGAGATGTTTTTCATCTATTATTTTAAAAACGTTGAAGGCGTTTTTAATTCTGCCTTCAACAGAGATCGACATACCTCTGATTATGATGTAA    ?A<AA<BDB?CDB@BCCCCCCDDDD@A@B<CCCEDCADDDDDDDCC@DDD@AABDDDDDDDDDDCEEDEEEFFFCBFFHGHHEH>CIJJJIIGJJGIIIGEIGD>HFGGDGGEGG@DIGGHFEFBCGHECGHHHIHHHHFDFEFFCC@    NM:i:5  MD:Z:50T31T36G4C18      MC:Z:30S118M    AS:i:116        XS:i:122        RG:Z:dad        XA:Z:chr20,+30424694,148M,6;chr4,-189273767,148M,8;chrUn_GL000219v1,-52002,12S136M,8;chr20,-29452832,12S136M,8;
```
The first column is the sequence name. The second is the "flag" column, which contains a variety of information about the alignment in a single number, all of which can be used for filtering reads. To decode this number, see [here](https://broadinstitute.github.io/picard/explain-flags.html). The third and fourth columns contain the reference sequence and position the read was mapped to. 

The fifth column is the [phred-scaled](https://en.wikipedia.org/wiki/Phred_quality_score) mapping quality (i.e. 10^( QUAL / -10) = probability of incorrect alignment). This is a measure of how certain the mapper is that the read belongs in this position in the reference. It is contingent on the reference genome, so if one region is missing from the reference, but similar to a second region that is present, reads derived from the first region may map with high quality to the second. 

After this comes the "cigar", which gives the actual alignment of the sequence to the reference in abbreviated format, followed by sequence and position of the mate pair sequence and the template length (if paired end sequencing was used). After that, the actual sequence and base qualities, followed by optional columns that may vary by read mapper. 


___

Now we can explore an alignment file. There are no scripts for the sections below. Run the code interactively. 

We can get a lot of basic stats on the SAM file using "samtools stats":

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools stats
samtools stats -r $GEN son.mkdup.bam >son.samstat.txt
```
If you do `less son.samstat.txt` you'll see this file is pretty messy, but we can pull out specific parts of it using `grep`. 

These are some basic stats about the alignment file:

```bash
grep ^SN son.samstat.txt | cut -f 2-
```

This is a histogram of per base coverage:

```bash
grep ^COV son.samstat.txt | cut -f 2-
```

And there is much more information. 

___

We can look at the per base coverage of individual regions easily using `samtools depth`

```bash
# make a list of bam files
ls *mkdup.bam >list.bam
# run samtools depth using that list
samtools depth -f list.bam -r chr20:29637000-29643000
```

We expect the parents to have 100x coverage and the son to have 50x coverage. Is there anything unusual about this region?

___

Samtools has implemented a simple alignment viewer in 'tview'. Let's look at the same region as above. 

```bash
# set a variable 'GEN' that gives the location of the reference genome:
GEN=/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK_all/resources_all/Homo_sapiens_assembly38.fasta
# run samtools tview
samtools tview --reference $GEN dad.mkdup.bam
```

After opening tview, type 'g' then enter 'chr20:29637000' to visit that location. 
"?" will bring up a help dialog
"q" exits

___

IGV, which has a graphical user interface, is a much better way of visualizing alignment files than `tview`. To use that, however, we'll have to run it locally. If you're working on the xanadu cluster, that means you'll have to download some files. 

For this exercise, we want to download `align_stepwise/son.mkdup.bam` and `align_stepwise/son.mkdup.bai`. You can do that using either the linux utility `scp` or using GUI ftp software. 

After you've done that and launched IGV, you'll need to change the default genome by going the "Genomes" menu, choosing "Load Genome From Server" and choosing "Human hg38". 

Then go to the "File" menu and choose "Load from File" and select son.mkdup.bam. 

You can visit the same region we examined previously, chr20:29637000-29643000, by entering the region in the dialog box at the top of the window and pressing the "Go" button. 

