#!/bin/bash

# run all tutorial scripts using slurm dependencies

jid1=$( sbatch --parsable 01_DataDownload/01_DownloadSamples.sh )
jid2=$( sbatch --parsable --dependency=afterok:$jid1  01_DataDownload/02_DownloadGenome.sh )

jid3=$( sbatch --parsable --dependency=afterok:$jid2  02_QualityControl/01_fastqc_raw.sh )
jid4=$( sbatch --parsable --dependency=afterok:$jid3  02_QualityControl/02_fastp.sh )
jid5=$( sbatch --parsable --dependency=afterok:$jid4  02_QualityControl/03_fastqc_post.sh )

jid6=$( sbatch --parsable --dependency=afterok:$jid5  03_AlignmentAndCoverage/01_bwa_index.sh )
jid7=$( sbatch --parsable --dependency=afterok:$jid6  03_AlignmentAndCoverage/02_bwa_align.sh )
jid8=$( sbatch --parsable --dependency=afterok:$jid7  03_AlignmentAndCoverage/03_index_bwa.sh )
jid9=$( sbatch --parsable --dependency=afterok:$jid8  03_AlignmentAndCoverage/04_samtools_index.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9  03_AlignmentAndCoverage/05_shortread_coverage.sh )
jid11=$(sbatch --parsable --dependency=afterok:$jid10 03_AlignmentAndCoverage/06_minimap2_align.sh )
jid12=$(sbatch --parsable --dependency=afterok:$jid11 03_AlignmentAndCoverage/07_longread_coverage.sh )
jid13=$(sbatch --parsable --dependency=afterok:$jid12 03_AlignmentAndCoverage/08_target_regions.sh )

jid14=$(sbatch --parsable --dependency=afterok:$jid13 04_VariantCalling_shortread/01_mpileup.sh )
jid15=$(sbatch --parsable --dependency=afterok:$jid14 04_VariantCalling_shortread/02_tabix.sh )
jid16=$(sbatch --parsable --dependency=afterok:$jid15 04_VariantCalling_shortread/03_freebayes.sh )
jid17=$(sbatch --parsable --dependency=afterok:$jid16 04_VariantCalling_shortread/04_gatk_gvcf.sh )
jid18=$(sbatch --parsable --dependency=afterok:$jid17 04_VariantCalling_shortread/05_gatk_DBimport.sh )
jid19=$(sbatch --parsable --dependency=afterok:$jid18 04_VariantCalling_shortread/06_gatk_genotype.sh )

jid20=$(sbatch --parsable --dependency=afterok:$jid19 05_VariantCalling_longread/01_clair3.sh )
jid21=$(sbatch --parsable --dependency=afterok:$jid20 05_VariantCalling_longread/02_convertcase.sh )
jid22=$(sbatch --parsable --dependency=afterok:$jid21 05_VariantCalling_longread/03_glnexus_clair3.sh )
jid23=$(sbatch --parsable --dependency=afterok:$jid22 05_VariantCalling_longread/04_gatk_single_sample_illumina.sh )

jid24=$(sbatch --parsable --dependency=afterok:$jid23 06_Filter_Compare/01_filter.sh )
jid25=$(sbatch --parsable --dependency=afterok:$jid24 06_Filter_Compare/02_summarize.sh )
jid26=$(sbatch --parsable --dependency=afterok:$jid25 06_Filter_Compare/03_compare.sh.sh )

jid27=$(sbatch --parsable --dependency=afterok:$jid26 07_Annotate/01_annotate_snpEff.sh )
jid28=$(sbatch --parsable --dependency=afterok:$jid27 07_Annotate/02_annotate_dbSNP.sh )