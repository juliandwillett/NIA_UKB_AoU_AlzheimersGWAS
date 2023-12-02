# While METAL is supposed to be able to handle different reference alleles, I was less confident of this (select variants were not captured)

awk 'BEGIN {OFS="\t"} NR==1 {print; next} { if ($6 > 0.5) { temp=$4; $4=$5; $5=temp; $9=-$9; $6=1-$6 } print }' \
/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_allminoralleles.txt

awk 'BEGIN {OFS="\t"} NR==1 {print; next} { if ($6 > 0.5) { temp=$4; $4=$5; $5=temp; $9=-$9; $6=1-$6 } print }' \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_allminoralleles.txt
