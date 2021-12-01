#!/bin/bash 
#SBATCH --job-name=mark_dups
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch/$USER

# son
IN=../results/align_stepwise/son.sort.bam
OUT=../results/align_stepwise/son.mkdup.bam

java -jar $PICARD MarkDuplicates \
        INPUT=$IN \
        OUTPUT=$OUT \
        METRICS_FILE=$IN.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=True

# mom
IN=../results/align_stepwise/mom.sort.bam
OUT=../results/align_stepwise/mom.mkdup.bam

java -jar $PICARD MarkDuplicates \
        INPUT=$IN \
        OUTPUT=$OUT \
        METRICS_FILE=$IN.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=True

# dad
IN=../results/align_stepwise/dad.sort.bam
OUT=../results/align_stepwise/dad.mkdup.bam

java -jar $PICARD MarkDuplicates \
        INPUT=$IN \
        OUTPUT=$OUT \
        METRICS_FILE=$IN.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        CREATE_INDEX=True

date