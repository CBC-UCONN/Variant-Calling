#!/bin/bash

# run all tutorial scripts using slurm dependencies

# run data download manually first, ensure the data were retrieved successfully, THEN run this script. 
# cd 01_DataDownload
# sbatch 01_DownloadSamples.sh
# sbatch 02_DownloadGenome.sh

# jid1=$( sbatch --parsable 01_DataDownload/01_DownloadSamples.sh )
# jid2=$( sbatch --parsable --dependency=afterok:$jid1  01_DataDownload/02_DownloadGenome.sh )

cd 02_QualityControl
jid3=$( sbatch --parsable   01_fastqc_raw.sh )
jid4=$( sbatch --parsable --dependency=afterok:$jid3  02_fastp.sh )
jid5=$( sbatch --parsable --dependency=afterok:$jid4  03_fastqc_post.sh )

cd ../03_AlignmentAndCoverage
jid6=$( sbatch --parsable --dependency=afterok:$jid5  01_bwa_index.sh )
jid7=$( sbatch --parsable --dependency=afterok:$jid6  02_bwa_align.sh )
jid8=$( sbatch --parsable --dependency=afterok:$jid7  03_index_bwa.sh )
jid9=$( sbatch --parsable --dependency=afterok:$jid8  04_samtools_index.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9  05_shortread_coverage.sh )
jid11=$(sbatch --parsable --dependency=afterok:$jid10 06_minimap2_align.sh )
jid12=$(sbatch --parsable --dependency=afterok:$jid11 07_longread_coverage.sh )
jid13=$(sbatch --parsable --dependency=afterok:$jid12 08_target_regions.sh )

cd ../04_VariantCalling_shortread
jid14=$(sbatch --parsable --dependency=afterok:$jid13 01_mpileup.sh )
jid15=$(sbatch --parsable --dependency=afterok:$jid14 02_tabix.sh )
jid16=$(sbatch --parsable --dependency=afterok:$jid15 03_freebayes.sh )
jid17=$(sbatch --parsable --dependency=afterok:$jid16 04_gatk_gvcf.sh )
jid18=$(sbatch --parsable --dependency=afterok:$jid17 05_gatk_DBimport.sh )
jid19=$(sbatch --parsable --dependency=afterok:$jid18 06_gatk_genotype.sh )

cd ../05_VariantCalling_longread
jid20=$(sbatch --parsable --dependency=afterok:$jid19 01_clair3.sh )
jid21=$(sbatch --parsable --dependency=afterok:$jid20 02_convertcase.sh )
jid22=$(sbatch --parsable --dependency=afterok:$jid21 03_glnexus_clair3.sh )
jid23=$(sbatch --parsable --dependency=afterok:$jid22 04_gatk_single_sample_illumina.sh )

cd ../06_Filter_Compare
jid24=$(sbatch --parsable --dependency=afterok:$jid23 01_filter.sh )
jid25=$(sbatch --parsable --dependency=afterok:$jid24 02_summarize.sh )
jid26=$(sbatch --parsable --dependency=afterok:$jid25 03_compare.sh.sh )

cd ../07_Annotate
jid27=$(sbatch --parsable --dependency=afterok:$jid26 01_annotate_snpEff.sh )
jid28=$(sbatch --parsable --dependency=afterok:$jid27 02_annotate_dbSNP.sh )