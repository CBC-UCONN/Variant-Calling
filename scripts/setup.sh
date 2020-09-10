#!/bin/bash

# run all tutorial scripts using slurm dependencies

jid1=$(sbatch --parsable Part1a_datadownload.sh )
jid2=$(sbatch --parsable --dependency=afterok:$jid1 Part1b1_fastqc.sh )
jid3=$(sbatch --parsable --dependency=afterok:$jid2 Part1b2_sickle_fastqc.sh )
jid4=$(sbatch --parsable --dependency=afterok:$jid3 Part1c_align.sh )
jid5=$(sbatch --parsable --dependency=afterok:$jid4 Part1d_compress.sh )
jid6=$(sbatch --parsable --dependency=afterok:$jid5 Part1e_sort.sh )
jid7=$(sbatch --parsable --dependency=afterok:$jid6 Part1f_markduplicates.sh )
jid8=$(sbatch --parsable --dependency=afterok:$jid7 Part1g_indexbams.sh )

jid9=$(sbatch --parsable --dependency=afterok:$jid8 Part2a_mpileup.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9 Part2b_variantcall.sh )
jid11=$(sbatch --parsable --dependency=afterok:$jid10 Part2c_tabix.sh )

jid12=$(sbatch --parsable --dependency=afterok:$jid2 Part3a_alignment.sh )
jid13=$(sbatch --parsable --dependency=afterok:$jid12 Part3b_indexbams.sh )

jid14=$(sbatch --parsable --dependency=afterok:$jid13 Part4a_coverage.sh )
jid15=$(sbatch --parsable --dependency=afterok:$jid14 Part4b_freebayes.sh )

jid16=$(sbatch --parsable --dependency=afterok:$jid15 Part4c_gatk_gvcf.sh )
jid17=$(sbatch --parsable --dependency=afterok:$jid16 Part4d_gatk_genomicsDBimport.sh )
jid18=$(sbatch --parsable --dependency=afterok:$jid17 Part4e_gatk_genotype.sh )

jid19=$(sbatch --parsable --dependency=afterok:$jid18 Part5a_filter.sh )
jid20=$(sbatch --parsable --dependency=afterok:$jid19 Part5b_compare.sh )

jid21=$(sbatch --parsable --dependency=afterok:$jid20 Part6a_annotate_SnpEff.sh )
jid22=$(sbatch --parsable --dependency=afterok:$jid21 Part6b_annotate_dbSNP.sh )
