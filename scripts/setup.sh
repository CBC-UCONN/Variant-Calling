#!/bin/bash

# run all tutorial scripts using slurm dependencies

# run data download manually first, ensure the data were retrieved successfully, THEN run this script. 
# cd 01_downloadData/
# sbatch 01_getSequences.sh
# sbatch 02_getGenome.sh

# jid1=$( sbatch --parsable 01_DataDownload/01_DownloadSamples.sh )
# jid2=$( sbatch --parsable --dependency=afterok:$jid1  01_DataDownload/02_DownloadGenome.sh )

cd 02_qc 
jid3=$( sbatch --parsable   01_fastqcRaw.sh )
jid4=$( sbatch --parsable --dependency=afterok:$jid3  02_trimmomatic.sh )
jid5=$( sbatch --parsable --dependency=afterok:$jid4  03_fastqcTrimmed.sh )

cd ../03_alignment
jid6=$( sbatch --parsable --dependency=afterok:$jid5  01_bwaIndex.sh )
jid7=$( sbatch --parsable --dependency=afterok:$jid6  02_bwaAlign.sh )

cd ../04_alignQC/
jid8=$( sbatch --parsable --dependency=afterok:$jid7  01_samstats.sh )
jid9=$( sbatch --parsable --dependency=afterok:$jid8  02_coverage.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9  03_bedtoolsNuc.sh )

cd ../05_variantCalling/
jid11=$(sbatch --parsable --dependency=afterok:$jid10 01_freebayes.sh )
jid12=$(sbatch --parsable --dependency=afterok:$jid11 02_bcftools.sh )
jid13=$(sbatch --parsable --dependency=afterok:$jid12 03_createDict.sh )
jid14=$(sbatch --parsable --dependency=afterok:$jid13 04_makeGVCFs.sh )
jid15=$(sbatch --parsable --dependency=afterok:$jid14 05_DBimport.sh )
jid16=$(sbatch --parsable --dependency=afterok:$jid15 06_genopeGVCFs.sh )

cd ../06_filteringAnnotating
jid17=$(sbatch --parsable --dependency=afterok:$jid16 01_filterVariants.sh )
jid18=$(sbatch --parsable --dependency=afterok:$jid17 02_normalizeVariants.sh )
jid19=$(sbatch --parsable --dependency=afterok:$jid18 03_dbSNP.sh )
jid20=$(sbatch --parsable --dependency=afterok:$jid19 04_bcftoolsCSQ.sh )
jid21=$(sbatch --parsable --dependency=afterok:$jid20 05_snpEff.sh )

# jid22=$(sbatch --parsable --dependency=afterok:$jid21 03_glnexus_clair3.sh )
# jid23=$(sbatch --parsable --dependency=afterok:$jid22 04_gatk_single_sample_illumina.sh )
# jid24=$(sbatch --parsable --dependency=afterok:$jid23 01_filter.sh )
# jid25=$(sbatch --parsable --dependency=afterok:$jid24 02_summarize.sh )
# jid26=$(sbatch --parsable --dependency=afterok:$jid25 03_compare.sh )
# jid27=$(sbatch --parsable --dependency=afterok:$jid26 01_annotate_snpEff.sh )
# jid28=$(sbatch --parsable --dependency=afterok:$jid27 02_annotate_dbSNP.sh )
