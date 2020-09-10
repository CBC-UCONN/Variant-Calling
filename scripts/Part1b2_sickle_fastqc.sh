#!/bin/bash
#SBATCH --job-name=trim_sickle
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=10G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname

module load sickle
module load fastqc

# trim son
SEQ=son
sickle pe -t sanger \
    -l 100 \
    -f ../rawdata/$SEQ.1.fq \
	-r ../rawdata/$SEQ.2.fq \
    -o ../rawdata/$SEQ.trim.1.fq \
    -p ../rawdata/$SEQ.trim.2.fq \
    -s ../rawdata/$SEQ.trim.0.fq
# trim mom
SEQ=mom
sickle pe -t sanger \
    -l 100 \
    -f ../rawdata/$SEQ.1.fq \
	-r ../rawdata/$SEQ.2.fq \
    -o ../rawdata/$SEQ.trim.1.fq \
    -p ../rawdata/$SEQ.trim.2.fq \
    -s ../rawdata/$SEQ.trim.0.fq
# trim dad
SEQ=dad
sickle pe -t sanger \
    -l 100 \
    -f ../rawdata/$SEQ.1.fq \
	-r ../rawdata/$SEQ.2.fq \
    -o ../rawdata/$SEQ.trim.1.fq \
    -p ../rawdata/$SEQ.trim.2.fq \
    -s ../rawdata/$SEQ.trim.0.fq


# run fastqc on all files matching "../rawdata/*trim*fq"
fastqc -t 6 -o ../fastqc/ ../rawdata/*trim*fq

